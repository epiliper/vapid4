to run: 


`docker build . -t <name_of_image>`

`docker run -it -v $(<filepath to dir with input files>) <name_of_image>`

`python3 -m vapid4 <input_fasta> <author .sbt file> --src_file <.src file> --output_loc <output dir> --align_dir <alignment dir>`
