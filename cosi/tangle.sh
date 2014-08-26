#!/usr/bin/env bash
# -*- mode: shell-script -*-
#
# tangle files with org-mode
#
DIR=`pwd`
FILES=""

# wrap each argument in the code required to call tangle on it
for i in $@; do
    FILES="$FILES \"$i\""
done

emacs -Q --batch \
--eval "(progn
(add-to-list 'load-path (expand-file-name \"/idi/sabeti-data/ilya/nsvn/Tools/org/devnew/org-mode/lisp/\"))
(add-to-list 'load-path (expand-file-name \"/idi/sabeti-data/ilya/nsvn/Tools/org/devnew/org-mode/contrib/lisp/\" t))
(require 'org)(require 'org-exp)(require 'ob)(require 'ob-tangle)
(mapc (lambda (file)
            (find-file (expand-file-name file \"$DIR\"))
            (org-babel-tangle)
            (kill-buffer)) '($FILES)))"

