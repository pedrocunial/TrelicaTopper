(TeX-add-style-hook
 "relatorio"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "paper=a4" "fontsize=11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("babel" "brazilian") ("microtype" "protrusion=true" "expansion=true") ("graphicx" "pdftex") ("appendix" "titletoc") ("adjustbox" "export")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "scrartcl"
    "scrartcl10"
    "inputenc"
    "fontenc"
    "fourier"
    "babel"
    "microtype"
    "amsmath"
    "amsfonts"
    "amsthm"
    "graphicx"
    "url"
    "indentfirst"
    "mathtools"
    "listings"
    "appendix"
    "float"
    "courier"
    "adjustbox"
    "sectsty"
    "fancyhdr")
   (TeX-add-symbols
    '("horrule" 1))
   (LaTeX-add-labels
    "sec:objetivos"
    "fig:trelica_entrada"
    "fig:secao14"
    "fig:ftool"
    "appendix:saida"
    "appendix:entrada")
   (LaTeX-add-bibliographies
    "referencias")
   (LaTeX-add-mathtools-DeclarePairedDelimiters
    '("Vector" "")))
 :latex)

