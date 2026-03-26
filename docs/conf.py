# Configuration file for the Sphinx documentation builder.

project = 'FEFF85 EXAFS C++'
copyright = '2025, Joint Photon Sciences Institute'
author = 'Joint Photon Sciences Institute'
release = '8.5.0'

extensions = [
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {
    'navigation_depth': 3,
    'collapse_navigation': False,
}
