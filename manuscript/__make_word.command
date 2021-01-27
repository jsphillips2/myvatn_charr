#!/bin/bash

## This creates the docx file from the __ms.tex document.
## On macOS, you should just be able to click on it from Finder

export name="charr_ms"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR"

pandoc ${name}.tex --bibliography=refs.bib -o ${name}.docx


exit
