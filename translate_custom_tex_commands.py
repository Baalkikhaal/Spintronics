import re


def match_pattern(pattern, string):

    # Pattern to match

    # pattern = "\\\\myScaSub{[^}]+}{[^}]+}"

    result = re.findall(pattern, string)
    
    return result

def iter_pattern(pattern, string):
    
    #pattern = "\\\\myScaSub{[^}]+}{[^}]+}"
    
    for match in re.finditer(pattern, string):
        print(match)
    
def substitute_pattern(pattern, replacement, string):
    #pattern = "\\\\myScaSub{[^}]+}{[^}]+}"
    
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

        # Return substitute
        substitute = replacement(match)
        substitutes.append(substitute)
    
    
    # Insert the translated substrings.
    
    joins = splits.pop(0)
    
    for split, substitute in zip(splits, substitutes):
        joins += substitute + split
    
    return joins

def replacement(match):
    """
    Return the substitute.
    """

    variables = re.findall("{[^}]+}", match)
    print(variables)

    substitute = "{}_{}".format(variables[0][1:-1], variables[1][1:-1])

    return substitute


def replacement1(match):
    """
    Return the substitute.
    """

    variables = re.findall("{[^}]+}", match)
    print(variables)

    substitute = "{}_{}".format(variables[0][1:-1], variables[1][1:-1])

    return substitute

def replacement2(match):
    """
    Return the substitute.
    """

    variables = re.findall("{[^}]+}", match)
    print(variables)

    #substitute = "\mathrm{}".format(variables[0][1:-1], variables[1][1:-1])
    # Use f-string to interpolate expressions.
    # Also to escape curly braces, use double braces.
    # To also evaluate expression within the braces, use triple braces.
    # Refer https://peps.python.org/pep-0498/#escape-sequences

    substitute = f' \\\\mathrm{{ \\\\mathbf{{ {variables[0][1:-1]} }} }} '
    return substitute

def replacement3(match):
    """
    Return the substitute.
    """

    variables = re.findall("[^\$]+", match)
    print(variables)

    substitute = f' \\\\( {variables[0]} \\\\) '

    return substitute

filename = "_unprocessed_notebooks/2022-04-08-Dispersion-in-refractive-index.ipynb"
#filename = "_notebooks/test.ipynb"

with open(filename, 'r') as f:
    stream = f.read()
    f.close()

#result = match_pattern(stream)

#print(result)

#iter_pattern(stream)

# list of patterns and corresponding replacement functions

#pattern = "\\\\myScaSub{[^}]+}{[^}]+}"

messages = [f'custom command \myScaSub',
            f'custom command \myVec',
            f'inline equation'
            ]
patterns = ["\\\\myScaSub{[^}]+}{[^}]+}",
            "\\\\myVec{[^}]+}",
            "\$[^\$\n]+\$"
            ]
replacements = [replacement1,
                replacement2,
                replacement3,
                ]

#joins = substitute_pattern(patterns[1], replacements[1], stream)

joins = stream

for message, pattern, replacement in zip(messages, patterns, replacements):
    print('='*40)
    print(f'finding and replacing {message}')
    print('='*40)
    joins = substitute_pattern(pattern, replacement, stream)
    stream = joins

fileout = "_notebooks/2022-04-08-Dispersion-in-refractive-index.ipynb"
#fileout = "_notebooks/test-pure-tex.ipynb"

with open(fileout, 'w') as f:
    f.write(joins)



