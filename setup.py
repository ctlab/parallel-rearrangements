from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='PaReBrick',
    version='0.1.0',
    description='A bioinf tool for finding genome rearrangements in bacterial genomes',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ctlab/parallel-rearrangements',
    author='Alexey Zabelkin',
    author_email='a.zabelkin@itmo.ru',
    classifiers=[
        'Topic :: Genome Rearrangements',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
    keywords='genome rearrangements, phylogenetic trees, non-convex characters, synteny blocks, phylogenetic tree, '
             'pattern consistency',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3.6, <4',
    install_requires=
        ['PyQt5', # for ete3 working properly
         'ete3', # phylogenetic trees
         'scikit-learn', # clustering characters
         'seaborn', # for drawer module
         'bg'], # breakpoint graphs
    entry_points={
        'console_scripts': [
            'PaReBrick=main:main',
        ],
    },
)