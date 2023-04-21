rule compute_answer:
    input:
        
    output:
        "src/tex/output/fit_table.txt"
    script:
        "src/scripts/fit.py"