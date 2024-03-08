to run: 

1. docker build . -t <name_of_image>

2. docker run -it -v $(<filepath to dir with input files>) <name_of_image>

3. python3 -m vapid4 <input_fasta> <author .sbt file> --src_file <.src file> --output_loc <output dir> --align_dir <alignment dir>
