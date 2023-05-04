import os

def write_trees_neightbours(blocks, ub_characters, neighbours, folder, show_branch_support, tree_holder, colors):
    for block, char in zip(blocks, ub_characters):
        labels = [f'{i}{"+" if i == len(colors) - 1 else ""} copies'
                  for i in range(max(char.values()) + 1)]

        file = folder + f'block_{block}.pdf'
        tree_holder.count_innovations_fitch(char)

        ns_block = neighbours[block]
        if len(set([frozenset(ns) for ns in ns_block.values()])) == 1: continue

        tree_holder.draw(file, legend_labels=labels, show_branch_support=show_branch_support, colors=colors,
                         neighbours=ns_block, neighbours_block=block, mode='r')