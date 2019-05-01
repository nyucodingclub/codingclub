List files in directory

ls **options** file

        	–l : list files with the permission, data, time, size

        	-h : human readable size

 

Sort files

sort **options** file

        	-n : numeric sort

        	-r : reverse sort

        	-k: position

        	-f : ignore case

 

Unique elements in file

uniq **options** file

        	-u : only lines that appear once

        	-d : only duplicated lines

        	-i : ignore case

        	-c :prefix line with number of occurences

 

Find pattern in file

grep **options** **pattern** file

        	-F : pattern is a list of fixed strings

        	-f : obtain pattern from file

        	-x : search for pattern that matches the whole line

        	-w : search for pattern that matches whole words

        	-i : ignore case

        	-v : find all lines NOT matching the pattern

        	-A NUM: get NUM lines after match

        	-B NUM : get NUM lines before match

        	-n : prefix each line of output with line number from input file

        	-E : interpret pattern as a regular expression

 

Replace pattern in file

sed **options** file

        	-i: replace within file (only works in linux..need to check)

        	sed 's/pattern/substitute/g' file  		 's'=substitute   'g'= globally

 

 

        	

 

 

