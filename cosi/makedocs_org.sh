#!/usr/bin/env bash
# -*- mode: shell-script -*-
#
# Generate cosi docs from Emacs Org-mode files.
#

emacs -Q --batch \
--eval "(progn
(add-to-list 'load-path (expand-file-name \"/idi/sabeti-data/ilya/nsvn/Tools/org/devnew/org-mode/lisp/\"))
(add-to-list 'load-path (expand-file-name \"/idi/sabeti-data/ilya/nsvn/Tools/org/devnew/org-mode/contrib/lisp/\" t))
(require 'org)(require 'org-exp)(require 'ob)(require 'ob-tangle)
(let ((org-publish-project-alist
       '(\"cosi\"
         :base-directory $srcdir :publishing-directory $builddir :publishing-function org-html-publish-to-html :include \"cosiguide.org\")))
   (org-html-publish-to-html)))"

