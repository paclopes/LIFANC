# Processes command line argument 1 to ease class use in matlab, namelly:
#
# 1 - Get the class constructor parameters from the first line of the class properties. 
# 2 - Add code to read the constructor parameters after the line with comment "% Read Parameters"
# 3 - Add code to copy the object properties to local variables after the line with comment "% Read properties".
# 4 - Add code to copy the object properties from local variables after the line with comment "% Write properties".
#
# Only updates the file when changes are detected.
#
# 22-6-2025 Paulo Lopes
# INESC-ID

import fileinput
import os
import sys

def update_file(file):
    properties_line = 0
    next_line = None

    for line in fileinput.input(file, inplace=True):
        if properties_line == 2:
            properties = line
            properties_line = 0;

            p = parameters.split(",") + properties.split(",")
            p = [x.strip() for x in p]

            properties_read_code = ""
            for i in range(0, len(p)):
                properties_read_code += p[i] + " = o." + p[i] + "; "

            properties_write_code = ""
            for i in range(0, len(p)):
                properties_write_code += "o." + p[i] + " = " + p[i] + "; "


        if properties_line == 1:
            parameters = line
            properties_line = 2

            p = parameters.split(",")
            p = [x.strip() for x in p]

            parameters_read_code = ""
            for i in range(0, len(p)):
                parameters_read_code += p[i] + " = param{" + str(i+1) + "}; "

        if not next_line is None:
            space = "";
            for char in line:
                if not char.isspace():
                    break
                space += char
            line = space + next_line;
            next_line = None

        if line.strip() == "properties":
            properties_line = 1

        if line.strip() == "% Read Parameters":
            next_line = parameters_read_code

        if line.strip() == "% Read properties":
            next_line = properties_read_code

        if line.strip() == "% Write properties":
            next_line = properties_write_code
        
        print(line.rstrip())

file = sys.argv[1]
time_file = file + "_time.txt"
if os.path.isfile(time_file):
    f = open(time_file)
    last_time = f.readline();
    f.close()
else:
    last_time = float('nan')
    
mtime = str(os.path.getmtime(file))

if last_time != mtime:
    update_file(file)
    f = open(time_file, "w")
    f.write(mtime)
    f.close()