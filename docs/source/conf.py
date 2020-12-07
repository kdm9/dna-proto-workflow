# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'NGS-Analysis-workflow v.1.0'
#copyright = '2020, PBGL'
author = 'Norman Warthmann (PBGL)'


# -- General configuration ---------------------------------------------------

# First need to install latex:
# sudo apt-get install  texmaker gummi texlive texlive-full texlive-latex-recommended latexdraw intltool-debian lacheck libgtksourceview2.0-0 libgtksourceview2.0-common lmodern luatex po-debconf tex-common texlive-binaries texlive-extra-utils texlive-latex-base texlive-latex-base-doc texlive-luatex texlive-xetex texlive-lang-cyrillic texlive-fonts-extra texlive-science texlive-latex-extra texlive-pstricks


# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
   'sphinx.ext.doctest',
   'sphinx.ext.intersphinx',
   'sphinx.ext.todo',
   'sphinx.ext.coverage',
   'sphinx.ext.mathjax',
   'sphinx.ext.ifconfig',
   'sphinx.ext.viewcode',
   'sphinx.ext.githubpages'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
source_suffix = '.rst'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

import sphinx_rtd_theme

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'default'
html_theme = 'sphinx_rtd_theme'
html_add_permalinks = ''
master_doc = 'index'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []


# -- Options for LaTeX output ---------------------------------------------

latex_engine = 'pdflatex'
latex_theme = 'howto'
latex_toplevel_sectioning = 'section'


# Configuration of Title Page
latex_maketitle = r'''
        \pagenumbering{Roman} %%% to avoid page 1 conflict with actual page

        \begin{titlepage}
            \vspace*{10mm} %%% * is used to give space from top
            \flushright\textbf{\Huge {PBGL-NGS-Analysis-Workflow v1.0}}

            \vspace{0mm} %%% * is used to give space from top
            \textbf{\Large {A Laboratory Manual}}

            \vspace{50mm}
            %%% * \textbf{\Large {Norman Warthmann}}

            \vspace{10mm}
            \textbf{\Large {Plant Breeding and Genetics Laboratory}}

            \vspace{0mm}
            \textbf{\Large {FAO/IAEA Joint Division}}

            \vspace{0mm}
            \textbf{\Large {Seibersdorf, Austria}}

	    \vspace{10mm}
            \normalsize Created: July, 2020

            \vspace*{0mm}
            \normalsize  Last updated: 24 November 2020

            %% \vfill adds at the bottom
            \vfill
            \small\flushleft {{\textbf {Please note:}} \textit {This is not an official IAEA publication but is made available as working material. The material has not undergone an official review by the IAEA. The views
expressed do not necessarily reflect those of the International Atomic Energy Agency or its Member States and remain the responsibility of the contributors. The use of particular designations of countries or territories does not imply any judgement by the publisher, the IAEA, as to the legal status of such countries or territories, of their authorities and institutions or of the delimitation of their boundaries. The mention of names of specific companies or products (whether or not indicated as registered) does not imply any intention to infringe proprietary rights, nor should it be construed as an endorsement or recommendation on the part of the IAEA.}}
        \end{titlepage}

        \pagenumbering{arabic}
        \newcommand{\sectionbreak}{\clearpage}
'''
latex_elements = {
   'releasename': 'Version 1.0',
   'maketitle': latex_maketitle,
}
