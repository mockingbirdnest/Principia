# Fold entries on one line and sort them alphabetically, thus by type then key.
([string]::join("`n", (get-content .\bibliography.bib)) `
     -replace "`n([^@])", ' $1').split("`n")            `
    | sort-object                                       `
    | out-file -encoding "UTF8" .\bibliography.bib
# Format with Biber.
biber --tool                         `
    --fixinits                       `
    --output-align                   `
    --output-indent=2                `
    --output-fieldcase=lower         `
    --output-file=.\bibliography.bib `
    .\bibliography.bib