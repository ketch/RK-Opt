"""Convert MATLAB code documentation to Restructured Text format suitable
   for use in Sphinx.

   This is used to generate the docs for the RK-opt package.

   Ignores additional functions in a single file.  Might be nice to
   change this behavior in the future.
"""
def extract_matlab_docstring(mfile):
    """Return the docstring from mfile, assuming that it consists of the
       first uninterrupted comment block.
    """
    docstring = ""
    instream = open(mfile,'rU')

    for i,line in enumerate(instream.readlines()[1:]):
        if line[0] == '%':
            if i>0:
                docstring += line[2:-1]+'\n'
            else:
                docstring += extract_function_name(line)+'\n'
                docstring += '='*len(line[2:])+'\n'
                docstring += '::\n\n'
                docstring += '    '+line[2:].rstrip()+'\n\n'
        else:
            return docstring+'\n'
    return docstring+'\n'

def compile_docstrings(directory,rstfile):
    """Write all the docstrings from directory to rstfile."""
    import os
    output=open(rstfile,'w')

    # Title
    write_h1(output,directory[3:-1])

    # Include the README, if there is one
    try:
        with open(os.path.relpath(directory+'README.rst')) as f:
            readme = f.read()
            output.write(readme)
            output.write('\n\n')
    except IOError:
        pass

    output.write('\n.. contents::\n\n')

    for fname in os.listdir(directory):
        if fname[-2:]=='.m':
            docstring = extract_matlab_docstring(os.path.relpath(directory+fname))
            output.write(docstring)
            output.write('\n\n')

def write_h1(stream,title):
    stream.write('.. _'+title+':\n\n')
    stream.write('='*len(title)+'\n')
    stream.write(title+'\n')
    stream.write('='*len(title)+'\n')

def write_h2(stream,title):
    stream.write(title+'\n')
    stream.write('='*len(title)+'\n')

def extract_function_name(signature):
    if '(' in signature and '=' in signature:
        i_start = signature.index('=')
        i_end = signature.index('(')
        return signature[i_start+1:i_end].lstrip()
    else:
        return signature.split()[-1]

if __name__ == '__main__':
    for subdir in ['RKtools','am_radius-opt','polyopt','RK-coeff-opt']:
        compile_docstrings('../'+subdir+'/',subdir+'.rst')
