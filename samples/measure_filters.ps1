Get-ChildItem -Path . -Filter *.fits | ForEach-Object {
    & python ..\measure_halfwidth.py -o "$($_.basename).csv" -w 40 -a 0.052 -c $_.basename $_.FullName
}