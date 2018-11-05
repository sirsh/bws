# Wecome

compile with *make* (see make file for simple compile line)

run using *bws.exe > out.data*

python *.\scripts\post_proc.py* can be used to process the out.data. It forces you to enter in case you accidentally overwrite files.
It is pretty basic starting point that assumes the output file is exactly called out.data and it asks for a name of an experiment to generate other files. 
One thing you might want to do is change it so it looks locally for out.data or change it to ask for the name of the data file to process.


# Params

There are various options as per help or see code in opts parsing. Note it is good to keep some sensible defaults in the code.

