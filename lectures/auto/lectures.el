;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "lectures"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("hyperref" "") ("graphicx" "") ("pgfpages" "") ("verbatim" "") ("bm" "") ("framed" "")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "hyperref"
    "graphicx"
    "pgfpages"
    "verbatim"
    "bm"
    "framed")
   (TeX-add-symbols
    "R"
    "jags"
    "wide"
    "bxt"
    "bx"
    "bxj"
    "bst"
    "bs"
    "bsi"
    "bsit"
    "bsitp"
    "ed"
    "cs"
    "dsixj"))
 :latex)

