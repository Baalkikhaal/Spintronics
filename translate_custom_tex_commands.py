import re


def match_pattern(string):

    # Pattern to match

    pattern = "\\\\myScaSub{[^}]+}{[^}]+}"

    result = re.findall(pattern, string)
    
    return result

def iter_pattern(string):
    
    pattern = "\\\\myScaSub{[^}]+}{[^}]+}"
    
    for match in re.finditer(pattern, string):
        print(match)
    
def substitute_pattern(string):
    pattern = "\\\\myScaSub{[^}]+}{[^}]+}"
    
    # Split the string at start of the matches
    splits = []
    m_starts = [match.start() for match in re.finditer(pattern, string)]
    m_ends = [match.end() for match in re.finditer(pattern, string)]
    
    start = 0
    
    for m_start, m_end in zip(m_starts, m_ends):
        split = string[start: m_start - 1]
        start = m_end
        splits.append(split)
    
    splits.append(string[start:])

    # Form the translated substring
    
    substitutes = []
    
    for match in re.findall(pattern, string):
        variables = re.findall("{[^}]+}", match)
        print(variables)
        substitute = "{}_{}".format(variables[0][1:-1], variables[1][1:-1])
        substitutes.append(substitute)
    
    
    # Insert the translated substrings.
    
    joins = splits.pop(0)
    
    for split, substitute in zip(splits, substitutes):
        joins += substitute + split
    
    return joins
    
filename = "_unprocessed_notebooks/2022-04-08-Dispersion-in-refractive-index.ipynb"
#filename = "_notebooks/test.ipynb"

with open(filename, 'r') as f:
    stream = f.read()
    f.close()

#result = match_pattern(stream)

#print(result)

#iter_pattern(stream)

joins = substitute_pattern(stream)

fileout = "_notebooks/2022-04-08-Dispersion-in-refractive-index.ipynb"
#fileout = "_notebooks/test-pure-tex.ipynb"

with open(fileout, 'w') as f:
    f.write(joins)



#translated = translate_pattern(stream, 

