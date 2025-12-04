# Encrypting data with Gocryptfs

Gocryptfs allows you to create a folder where the file contents and names are encrypted before being saved to the underlying storage system.  This can be used to allow you to store data in places you might otherwise not be able to.

Technical details on Gocryptfs can be found on their site at [https://github.com/rfjakob/gocryptfs](https://github.com/rfjakob/gocryptfs)

# Temporarily store data
This example explains how you can temporarily download and use data that you otherwise could not have in the HPC environment.  

## Process Description
* Allocate a HPC compute node using Slurm and connect to it
* Setup a temporary encrypted folder within temp space on the compute node.
* Download your sensitive data from its original location, in this example Box, directly into this encrypted folder.
* Run your analysis
* Copy your completed analysis results to another location.  (Make sure that location is approved for the data use)
* Delete your temporary encrypted folder 
* Disconnect from compute node

## Technical steps

This example creates a folder in tmp to store encrypted contents using a temporary key and password.  The key, password and data is deleted upon job completion.  

```
module load gocryptfs
mkdir ~/securetemp
cd ~/securetemp

temp_encrypted_dir=$(mktemp -d)

mkdir access 
temp_password_file=$(mktemp)
head /dev/urandom | tr -dc A-Za-z0-9_ | head -c 32 > ${temp_password_file}
chmod 600 temp_password_file
#Make  it
~/gocryptfs/gocryptfs -init ${temp_encrypted_dir} -i 60m -passfile ${temp_password_file}
#Mount it
~/gocryptfs/gocryptfs ${temp_encrypted_dir} access -passfile ${temp_password_file}
rm ${temp_password_file}  #After this is run the dir can not be mounted any additional times.  
cd ~/securetemp/access

#DO YOUR SECURE WORK
#Making sure to only write confidential data to ~/securetemp/access
#Remember to copy results somewhere?  Anything in here will be deleted.


#AFTER WORK COMPLETE
#Remove all files when complete
rm -rf ~/securetemp
rm -rf ${temp_encrypted_dir}
```

