BaRC code for a variety of bioinformatics tasks, written by members of the 
[Whitehead Institute](https://wi.mit.edu/) [Bioinformatics and Research Computing](http://barc.wi.mit.edu/) group.

For all code, executing the script by itself should return a summary and a usage statement.
Each directory includes code, sample input and output, and sample commands (in README).

Code is organized by programming language, but since all scripts are designed to be run on the command line, independently of language, this division is not so practically meaningful.

To get access to all of the BaRC code, download the code as a ZIP archive, and then one way to proceed is

```
# Unzip
unzip barc-main.zip
# Enter directory
cd barc-main
# Make all code executable
find Perl   -name \*.pl -exec chmod 755 {} \;
find Python -name \*.py -exec chmod 755 {} \;
find R      -name \*.R -exec chmod 755 {} \;
find shell  -name \*.sh -exec chmod 755 {} \;
# Create a new directory for links to all executables
mkdir BaRC_code
# Add links to all executables
find $(pwd)/Perl   -name \*.pl -executable -exec ln -sf {} BaRC_code \;
find $(pwd)/Python -name \*.py -executable -exec ln -sf {} BaRC_code \;
find $(pwd)/R      -name \*.R  -executable -exec ln -sf {} BaRC_code \;
find $(pwd)/shell  -name \*.sh -executable -exec ln -sf {} BaRC_code \;
# Add BaRC_code (with executables) to PATH
export PATH=$PATH:$(pwd)/BaRC_code
```
