#!/usr/bin/env Rscript

#biasv=c(0.5,0.6,0.7,0.8,0.9,0.99)
biasv = c(0.7)
nrep =20 

maxgen = 10000

a <- c(0.8)

b <- c(0.0025)

c <- c(0)

# generate a date_time stamp as a character
date_str <- format(Sys.time(), "%d%m%Y_%H%M%S")

# generate names for the output file based on date and time 
# to ensure some sort of uniqueness, to avoiding that output
# files from different runs overwrite eachother
output_file_prefix = paste0("2sim_good_genes_",date_str)

counter <- 0

exe = "./good_genes.exe"

batch_file_contents <- ""

for (rep_i in 1:nrep)
{
    for (biasv_i in biasv)
    {
	    for (a_i in a)
	    {
		    for (b_i in b)
		    {
			   for (c_i in c)

        counter <- counter + 1
        file_name_i <- paste0(output_file_prefix,"_",counter)

        echo_str <- paste("echo",counter)

        command_str <- paste(exe,
                        biasv_i,
			a_i,
			b_i,
			c_i,
                        format(maxgen,scientific=F),
                        file_name_i)

        # append to batch file contents
        batch_file_contents <- paste0(batch_file_contents
                ,"\n"
                ,echo_str
                ,"\n"
                ,command_str)
		    }
	    }
    }
}

writeLines(text=batch_file_contents)
