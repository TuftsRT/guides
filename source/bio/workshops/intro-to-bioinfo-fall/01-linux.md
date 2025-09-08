# What is Linux?

Linux is a free, open-source, and Unix-like operating system kernel that was originally created by **Linus Torvalds** in 1991. Over time, Linux has grown into a full-fledged operating system used worldwide across various types of devices, from servers and desktop computers to smartphones and embedded systems.

**Popular Linux Distributions:**

​ • **Ubuntu:** A user-friendly distribution popular for desktop and server use, based on Debian.

​ • **Fedora:** A cutting-edge distribution often used by developers and those who want the latest features.

​ • **Debian:** Known for its stability and extensive software repositories, often used in server environments.

​ • **CentOS/AlmaLinux/Rocky Linux:** Enterprise-grade distributions derived from Red Hat Enterprise Linux (RHEL).

​ • **Arch Linux:** A rolling release distribution known for its simplicity and customization, aimed at advanced users.

​ • **Kali Linux:** A distribution designed for penetration testing and security research.

<img src="https://nixwindows.wordpress.com/wp-content/uploads/2015/02/linux-distro-stickers.png" width="600">


# Our clusters' OS
![](img/os.png)

## Files and File System

### Everything is a file

**A file is an addressable location that contains some data which can take many forms.**

- Text data

- Binary/Image data

**Files have associated meta-data**

- Owner

- Group

- Timestamps

- Permission:

  1. Read r
  1. Write w
  1. Execute x
  1. no permission -

### File permissions

<img src="https://i.imgur.com/yxNrpKJ.png" width="600">

### File organization

Everything is mounted to the root directory

Files are referred to by their location called **path**

- Absolute Path (From the root): **/cluster/tufts/mylab/user01**

- Relative Path (From my current location)：**user01**

## Must-known Linux/Unix Tools

### File and Directory Management

Linux provides powerful tools for managing files and file systems. Here we will introduce a few essential commands.

#### pwd: print the current working directory

##### Usage

```
$ pwd
/cluster/home/yzhang85
$ cd /cluster/tufts/rt/yzhang85/
$ pwd
/cluster/tufts/rt/yzhang85
```

#### cd: change directory

##### Usage

```
cd [directory]
```

If a directory is not supplied as an argument, it will default to your **home** directory.

```
$ pwd
/cluster/tufts/rt/yzhang85
$ cd ..
$ pwd
/cluster/tufts/rt
$ cd
$ pwd
/cluster/home/yzhang85
```

##### Shortcuts

- **..**: cd to the parent directory. **cd ..**
- **~**: cd to the home directory. **cd ~**
- **-**: cd to the previous directory. **cd -**

#### ls: list all the files in the given directory

##### Usage

```
ls [options] [directory]
```

##### Common options:

- **-1**: list each file/directory on a separate line
- **-l**: lists files/directories with their most common metadata
- **-a**: include hidden files /directories (files’ name begins with a dot **.**)
- **-h**: print size of files/directories in human readable format

#### chmod: manage file permissions

##### Symbolic Notation

- **u**: user (owner)
- **g**: group
- **o**: others
- **a**: all (user, group, others)
- **+**: add permission
- **-**: remove permission
- **=**: set exact permission

##### Examples

```
$ chmod g+w filename ## Give the group write permission
$ chmod u+x filename ## Give user execute permission
$ chmod a+r filename ## Give all users read access
$ chmod u=rw,g=r,o=r filename ## Give user read and write permission, group and other only read permission.
```

##### Recursive updating permissions with -R

To apply permissions recursively to all files and subdirectories within a directory, use the **-R** option:

```
$ chmod -R g+rx /path/to/directory
```

#### touch: create new files or update timestamps

touch is used to create new files or to update the timestamps (access and modification times) of existing files.

##### Create new file

```
touch newfile.txt
```

##### Update timestamps of existing files

```
touch existingfile.txt
```

#### mkdir: create new directory

##### Usage

```
mkdir [options] dir_name
```

##### Common option

- **-p**: Creates parent directories if they don't exist.

```
$ mkdir -p rnaseq/output
```

This will create `output` folder as well as its parent folder `rnaseq` if it doesn't exist.

#### mv: move a file/directory to a new location or rename it

##### Usage

```
mv [options] source destination
```

##### Common option

- **-i**: Prompts for confirmation before overwriting an existing file. Useful to avoid accidental data loss.
- **-f**: Forces the operation without prompting, even if an existing file would be overwritten. **Use with caution**!

#### cp: copy a file/directory

##### Usage

```
cp [options] source destination
```

##### Common option

- **-r**: To copy directory

#### rm: remove files/directories

##### Usage

```
rm [options] file/directory
```

##### Common option

- **-r:** Deletes recursively any file and subdirectories contained within the given directory.

### Text processing

Linux command-line tools are invaluable for bioinformatics text processing due to their efficiency and flexibility. They allow for rapid manipulation and analysis of large biological datasets, such as DNA sequences, protein structures, and gene expression data. Commands like `grep`, `sed`, `awk`, and `cut` are essential for filtering, extracting, and reformatting text-based biological information.

#### cat: catenate files(joins their contents)

##### Usage

```
cat [options] file1 file2 …
```

##### Common option

- **-n:** tells cat to **number each line of the output**. This is helpful for debugging scripts.

#### head/tail: display the beginning/end of a file

##### Usage

```
head/tail [options] file
```

##### Common option

- **-n** \[number\]: Specifies the number of lines to display (default: 10).

#### less/more: view the content of a file page by page

##### Usage

```
less largefile.txt
more largefile.txt
```

##### Navigation with less

- **Arrow Keys**: Use the up and down arrow keys to scroll line by line.
- **Spacebar**: Move forward one page.
- **b**: Move backward one page.
- **q**: Quit and exit less.
- **/search_term**: Search for search_term within the file. Press n to go to the next occurrence.
- **g**: Go to the beginning of the file.
- **G**: Go to the end of the file.

##### Navigation with more

- **Spacebar**: Move forward one page.
- **Enter**: Move forward one line.
- **q**: Quit and exit more.
- **/search_term**: Search for search_term within the file (forward only).

#### cut: Extract sections from each line of files

cut is a text-processing utility used to extract specific sections of lines from a file or standard input. It is commonly used to split lines of text based on delimiters, extract columns of data, and work with fixed-width fields. The cut command is helpful for manipulating text files like **TSVs** and **CSVs**.

##### Usage

```
cut [OPTIONS] [FILE...]
```

##### Common options

- **-d**: Specifies the delimiter for field extraction. Default is TAB. Example: **-d,** (use a comma as the delimiter).
- **-f**: Selects specific fields, used with a delimiter (default is TAB). Example: **-f 1,3** (extract fields 1 and 3).

##### Example

```
cut -f1,3 -d, file.csv ##(Extract columns 1 and 3 from a comma-separated file)
```

#### sort: Sort lines of text files

`sort` is designed to sort plain-text data with columns. Running **sort** without any arguments simply sorts a file alphabetically.
By default, sort treats **blank characters** (like tab or spaces) as field delimiters. If your file uses another delimiter (such as a comma for CSV files), you can specify the field separator with **-t** (e.g., **-t ","**).

##### Usage

###### Sort lines alphabetically

```
sort file.txt
```

###### Sort lines numerically

```
sort -n file.txt
```

###### Sort by a specific column

```
sort -k 2 file.txt  # Sort by the second column
```

###### Sort by multiple columns

```
sort -k 1,2 file.txt  # Sort by the first column, then the second
```

#### uniq: Report or filter out repeated lines in a file.

uniq is used to filter out or identify duplicate lines in a sorted file. It removes adjacent duplicate lines, so it’s often combined with the `sort` command to find unique lines in a file effectively.

##### Usage

```
sort file.txt | uniq
```

#### grep:Extracting lines matching (not matching) a pattern

##### Usage

```
grep [options] PATTERN file
```

##### Common option

- **-i**: ignore cases
- **-v**: select non-matching lines.
- **-A NUM:** Print **NUM** lines of trailing context after matching lines.
- **-B NUM:** Print **NUM** lines of leading context before matching lines.

#### sed: Stream editor for modifying file content

sed (short for stream editor) is a powerful text-processing tool in Bash that allows you to parse and transform text in files or streams. It is commonly used to perform basic text manipulations like search and replace, insert and delete lines, and apply regular expressions on text data.

##### Substitution (Search and Replace)

Replace the first occurrence of **old** with **new** in each line:

```
sed 's/old/new/' filename.txt
```

Replace **all** occurrences of **old** with **new** in each line:

```
sed 's/old/new/g' filename.txt
```

##### In-place substitution

```
sed -i 's/old/new/g' filename.txt
```

Warning: Use this command with caution as it directly modifies the original file. To create a backup, use `-i.bak`:

```
sed -i.bak 's/old/new/g' filename.txt
```

##### Delete lines

```
sed '/pattern/d' filename.txt
```

`sed` and `awk` that we will introduce later are very powerful bash commands. If you have interest in learning more about their usage, below is a very good book.
<img src="https://m.media-amazon.com/images/I/916OGXYyYVL._SL1500_.jpg" width="400">

### Data Compression and Archiving

When working with files on Linux, compressing them to save space and bundling multiple files into a single archive is a common practice. The commands gzip, gunzip, and tar are essential tools for file compression and archiving in Bash.

#### gzip: Compress files

gzip is used to compress files in Linux. It reduces file sizes using the **DEFLATE** algorithm, resulting in files with the `.gz` extension.

##### Usage

###### Compress file/files

```
gzip file.txt ## Compress file.txt and replace it with file.txt.gz.
gzip file1.txt file2.txt file3.txt ## Each file will be compressed and replaced with a .gz version.
```

###### Keep the Original File

```
gzip -k file.txt ## Create file.txt.gz while preserving the original file.txt.
```

###### Compress the folder

```
gzip -r directory_name
```

#### gunzip: decompress .gz files

##### Decompress files

```
gunzip file.txt.gz
gunzip file1.txt.gz file2.txt.gz
```

##### Keep the Compressed File

```
gunzip -k file.txt.gz
```

#### tar: Archive multiple files into one or extract them.

tar is used to create, extract, and manipulate archive files. However, tar itself does not compress files; it only archives them by combining multiple files and directories into a single file. This file usually has a `.tar` extension. However, tar can be used in combination with other compression utilities (like `gzip` or `bzip2`) to compress the archive.

##### Create a .tar archive without compression

```
tar -cvf archive.tar my_folder
```

##### Extract a tar file

```
tar -xvf archive.tar
```

##### Creating a compressed archive(.tar.gz)

```
tar -cvzf archive.tar.gz my_folder
```

##### Extracting a compressed archive(.tar.gz)

```
tar -xvzf archive.tar.gz
```

### Other useful tools

#### Environment variables

##### Define variables

```
VARIABLE=value    ## No space around =
```

##### Variable reference

```
$VARIABLE     ## echo $VARIABLE
```

##### Commonly used environment variables

- **\$USER**: the login user
- **\$HOME**: the home directory
- **\$PWD**: the current directory
- **\$PATH**: A list of directories which will be checked for executable files

#### Redirection: >, >>, \<

- `>`: Overwrites the contents of a file with the command's output

  ```
  `cat file1 file2 > files`
  ```

- `>>`: Appends the output to the end of an existing file

​ `cat file3 >> files`

- `<`: Uses the contents of a file as input to a command

​ `sort < names.txt`

#### Pipe: |

Pipes in Linux are a powerful feature that allows you to connect the output of one command directly as the input to another command. This is a key concept in Unix/Linux philosophy, which promotes the use of small, modular tools that can be combined to perform complex tasks.

A pipe is represented by the `|` symbol. When you place a pipe between two commands, the standard output (`stdout`) of the command on the left of the pipe becomes the standard input (`stdin`) for the command on the right.

##### Usage

```
command1 | command2
```

```
sort file.txt | uniq
```

- **sort file.txt**: Sorts the lines in file.txt.

- **uniq**: Removes duplicate lines from the sorted output.

#### Wildcards: selecting multiple files/directories based on patterns

- **\***: Represents zero or more characters:

  ```
   - **\*.fastq.gz**  matches all fastq.gz files
  ```

- **?**: Represents a single character:

  ```
   - **file?.txt** matches "file1.txt", "fileA.txt", but not "file12.txt".
  ```

- **[]**: Represents a single character within a specified range or set:

  ```
   - **[abc]at** matches "bat", "cat", or "aat".
   - **[0-9]** matches any single digit.
  ```

#### Alias

An alias in Linux is a custom shortcut or abbreviation for a command or a series of commands. Once defined, you can use the alias in place of the original command.

##### Creating an alias

To create an alias, use the alias command followed by the name you want to give the alias and the command it should execute.

```
alias alias_name='command'
```

##### Example

```
alias ll='ls -l'
alias la='ls -a'
alias mav='module avail'
alias ml='module load'
```

##### Using alias

```
[yzhang85@login-prod-03 MPI]$ ll
total 31
-rwxrwx--- 1 yzhang85 yzhang85 16256 May 16 07:57 hello_mpi*
-rw-rw---- 1 yzhang85 yzhang85   731 May 16 07:56 hello_mpi.c
-rwxrwxr-x 1 yzhang85 yzhang85  8392 Dec 23  2023 hello_world*
-rwxrwxr-x 1 yzhang85 yzhang85   731 Oct 31  2023 hello_world.c*
-rwxrwxr-x 1 yzhang85 yzhang85  4053 Oct 31  2023 mpi_hello_world.c*

[yzhang85@login-prod-03 MPI]$ mav blast
-------------------------------- /opt/shared/Modules/modulefiles-rhel6 --------------------------
   blast/2.2.24    blast/2.2.31    blast/2.3.0    blast/2.8.1

------------------------------- /cluster/tufts/hpc/tools/module ---------------------------------
   blast-plus/2.11.0    ncbi-magicblast/1.5.0

------------------------------- /cluster/tufts/biocontainers/modules -----------------------------
   blast/2.15.0 (D)

  Where:
   D:  Default Module
```

#### ln -s: softlink

A **soft link** (also known as a symbolic link or symlink) is a type of file in Linux that points to another file or directory. It’s essentially a shortcut that references the location of another file, allowing you to access it from a different location in the filesystem.

##### Usage

To create a soft link, you use the `ln` command with the `-s` option:

```
ln -s target_file link_name
```

- target_file: The file or directory you want to link to.
- link_name: The name of the symlink that will point to the target.

## Conditional Statements and Loops

### for

Using `for` loop, we can execute a set of commands for a finite number of times for every item in a list.

#### Syntax

```
for VARIABLE in LIST
do
  command1
  command2
done
```

#### Example

```
for query in /cluster/tufts/Yourlab/input/*.fasta
do
    echo "Running BLAST for $query"
    blastn -query $query -db /path/to/database -out ${query%.fasta}_blast_results.txt -outfmt 6
done
```

### while

The command next to `while` is evaluated. If it is evaluated to be true, then the commands between `do` and `done` are executed.

#### Syntax

```
while [ condition ]
do
   command1
   command2
   command3
done
```

#### Example

```
#!/bin/bash

# File containing the list of sequence files
file_list="sequence_files.txt"
# Path to the BLAST database
blast_db="/path/to/blast_database"

# Read the list of files line by line
while IFS= read -r sequence_file
do
    # Check if the file exists
    if [ -f "$sequence_file" ]; then
        echo "Processing $sequence_file"
        # Run BLAST on the current sequence file
        blastn -query "$sequence_file" -db "$blast_db" -out "${sequence_file%.fasta}_blast_results.txt" -outfmt 6
        echo "Finished processing $sequence_file"
    else
        echo "File $sequence_file not found"
    fi
done < "$file_list"

echo "All files have been processed."
```

The content of `sequence_files.txt` is the name of all of the query fasta files:

```
sequence1.fasta
sequence2.fasta
sequence3.fasta
```

### if: conditionals in a bash script

In bioinformatics, if statements are used to make conditional decisions in scripts, such as checking whether a file or result already exists before performing data analysis or processing tasks, thus optimizing workflows and avoiding redundant computations.

#### Syntax

```
if [ commands ]
then
  [ if-statements ]
elif [commands]
  [ elif-statements ]
else
  [ else-statements ]
fi
```

#### Example

```
#!/bin/bash

# Define the paths for input files and output file
input_bam="/path/to/genomic_data.bam"
output_vcf="/path/to/output_directory/variants.vcf"

# Check if the output VCF file already exists
if [ -f "$output_vcf" ]; then
    echo "Output file $output_vcf already exists. Skipping the variant calling step."
else
    echo "Output file $output_vcf does not exist. Running variant calling."

    # Run the variant calling tool (e.g., bcftools)
    bcftools mpileup -Ou -f /path/to/reference_genome.fasta "$input_bam" | \
    bcftools call -mv -Ov -o "$output_vcf"

    echo "Variant calling complete. Results saved to $output_vcf"
fi
```

This if statement ensures that the variant calling process is skipped if the output file already exists, preventing unnecessary computations and saving processing time. This approach is useful in bioinformatics workflows where redundant analyses can be avoided by checking for existing results.
