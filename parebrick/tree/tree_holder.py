from ete3 import Tree, TreeStyle, TextFace, RectFace, SeqMotifFace

from collections import defaultdict, Counter
from itertools import combinations


def get_neighbour_motifs(n, ns_colors, offset, inverse=False):
    block, orient = n[:-1], n[-1]
    clr = ns_colors[block]

    if (orient == 'h' and not inverse) or (inverse and orient == 't'):
        return [[offset, offset + 30, "[]", None, 10, clr, clr, f"arial|2|white|+{block}"],
                [offset + 30, offset + 40, ">", None, 10, clr, clr, None]]
    else:
        return [[offset, offset + 10, "<", None, 10, clr, clr, None],
                [offset + 10, offset + 40, "[]", None, 10, clr, clr, f"arial|2|white|-{block}"]]


def generate_neighbour_face(node_ns, ns_colors, block):
    motifs = [[0, 0, "blank", None, 10, None, None, None]]

    for i, (nbr1, nbr2) in enumerate(sorted(node_ns)):
        offset = 160 * i

        motifs.extend(get_neighbour_motifs(nbr1, ns_colors, offset))
        motifs.extend(get_neighbour_motifs(f'{block}h', ns_colors, offset + 50))
        motifs.extend(get_neighbour_motifs(nbr2, ns_colors, offset + 100, True))

        motifs.append([offset + 140, offset + 160, "blank", None, 10, None, None, None])

    return SeqMotifFace('', motifs=motifs)


class TreeHolder:
    def __init__(self, tree, logger, scale=None, labels_dict=None, node_colors=defaultdict(lambda: 'black')):
        self.tree = Tree(tree)
        self.scale = scale

        for node in self.tree.traverse():
            if len(node.children) == 3:
                logger.info("Trying to root tree by first child of root")
                logger.info(f'Children of root: {node.children}')
                self.tree.set_outgroup(node.children[0])
            break

        for node in self.tree.traverse():
            # Hide node circles
            node.img_style['size'] = 0

            if node.is_leaf():
                try:
                    name_face = TextFace(labels_dict[node.name] if labels_dict else node.name,
                                         fgcolor=node_colors[node.name])
                except KeyError:
                    msg = f'There is not label for leaf {node.name} in labels file'
                    logger.error(msg)
                    raise KeyError(msg)
                node.add_face(name_face, column=0)

    def draw_neighbours(self, neighbours, block, colors=['Crimson', 'Teal', 'DarkGreen', 'Purple', 'DarkKhaki',
                                                         'MediumVioletRed', 'DarkOrange', 'Navy', 'RosyBrown',
                                                         'DarkGoldenrod', 'Sienna', 'Indigo', 'DarkRed', 'Olive', 'Black']):
        posible_ns = sorted(list(set(n[:-1] for nss in neighbours.values() for ns in nss for n in ns)))

        ns_colors = {posible_ns[i]: colors[i % len(colors)] for i in range(len(posible_ns))}
        ns_colors[str(block)] = 'grey'

        # if block != 2: return

        for node in self.tree.traverse():
            if not node.is_leaf(): continue
            face = generate_neighbour_face(sorted(neighbours[node.name]), ns_colors, block)
            node.add_face(face, 1, "aligned")

    def draw(self, file, colors, color_internal_nodes=True, legend_labels=(), show_branch_support=True,
             show_scale=True, legend_scale=1, mode="c", neighbours=None, neighbours_block=None):
        max_color = len(colors)

        used_colors = set()
        for node in self.tree.traverse():
            if not (color_internal_nodes or node.is_leaf()): continue
            color = colors[min(node.color, max_color - 1)]
            node.img_style['bgcolor'] = color
            used_colors.add(color)

        ts = TreeStyle()
        ts.mode = mode
        ts.scale = self.scale
        # Disable the default tip names config
        ts.show_leaf_name = False
        ts.show_branch_support = show_branch_support

        # ts.branch_vertical_margin = 20
        ts.show_scale = show_scale
        cur_max_color = max(v.color for v in self.tree.traverse())
        current_colors = colors[0:cur_max_color + 1]

        for i, (label, color_) in enumerate(zip(legend_labels, current_colors)):
            if color_ not in used_colors: continue
            rf = RectFace(20 * legend_scale, 16 * legend_scale, color_, color_)
            rf.inner_border.width = 1
            rf.margin_right = 14
            rf.margin_left = 14

            tf = TextFace(label, fsize=26 * legend_scale)
            tf.margin_right = 14

            ts.legend.add_face(rf, column=0)
            ts.legend.add_face(tf, column=1)

        if neighbours:
            old_tree = self.tree.copy()
            self.draw_neighbours(neighbours, neighbours_block)

        self.tree.render(file, w=1000, tree_style=ts)

        if neighbours:
            self.tree = old_tree

    def get_all_leafs(self):
        return {node.name for node in self.tree.get_leaves()}

    def count_innovations_fitch(self, leaf_colors, count_second_color=True):
        def assign_colorset_feature(v):
            if v.is_leaf():
                v.add_features(colorset={leaf_colors[v.name]}, color=leaf_colors[v.name])
            else:
                try:
                    child1, child2 = v.children
                except ValueError:
                    print(v.children)
                    raise ValueError('Tree must me binary')
                cs1 = assign_colorset_feature(child1)
                cs2 = assign_colorset_feature(child2)
                v.add_features(colorset=(cs1 & cs2) if len(cs1 & cs2) > 0 else cs1 | cs2)

            return v.colorset

        def chose_color(colorset):
            return sorted(colorset, key=lambda c: color_counter[c], reverse=True)[0]

        def down_to_leaves(v, color):
            if v.is_leaf(): return
            v.add_features(color=color if color in v.colorset else chose_color(v.colorset))
            for child in v.children:
                down_to_leaves(child, v.color)

        def count_innovations(v, innovations):
            for child in v.children:
                if v.color != child.color and not (not count_second_color and (v.color == 2) or (child.color == 2)):
                    innovations[child.color].append(child)
                count_innovations(child, innovations)

        color_counter = Counter(leaf_colors.values())

        # get colorsets for internal nodes
        root = self.tree.get_tree_root()
        assign_colorset_feature(root)

        # get color for internal nodes
        root_color = chose_color(root.colorset)
        down_to_leaves(root, root_color)

        # get inconsistent colors
        self.innovations = defaultdict(list)
        count_innovations(root, self.innovations)

    def count_parallel_rearrangements(self, skip_grey):
        score, count, count_all = 0, 0, 0
        for color, nodes in self.innovations.items():
            if len(nodes) <= 1 or (skip_grey and color == 1): continue
            count += 1
            count_all += len(nodes)
            for n1, n2 in combinations(nodes, 2):
                score += n1.get_distance(n2)
        return score, count, count_all

    def count_parallel_breakpoints(self):
        count = sum(map(len, self.innovations.values()))
        score = sum(
            n1.get_distance(n2) for n1, n2 in combinations((n for ns in self.innovations.values() for n in ns), 2))
        return score, count

    def draw_coloring(self, file):
        for node in self.tree.traverse():
            node.img_style['bgcolor'] = self.colors[node.color]
        ts = TreeStyle()
        ts.show_leaf_name = False
        self.tree.render(file, w=1000, tree_style=ts)

    def prune(self, ls):
        self.tree.prune(list(ls))
