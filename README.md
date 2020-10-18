# PhyloPy
A data driven protein ortholog finding tool utilizing phylogenetic trees, intron/exon boundaries, genomic context, and protein domain figures

# Setting Up Conda Environment

  1. Click [here](https://www.anaconda.com/products/individual) and install Anaconda for your appropriate system
  
  
  2. Navigate to your downloads folder and run the shell script
    
      For python 3.7: `bash ~/Downloads/Anaconda3-2020.02-Linux-x86_64.sh`

      For python 2.7: `bash ~/Downloads/Anaconda2-2019.10-Linux-x86_64.sh`
      
      
  3. Create the virtual environment in the cloned directory (see Cloning the Github Repo). Make sure evironment.yml file is in dir.
  
      `conda PhyloPy create -f environment.yml python=2.7`
      
  4. Install the python dependencies
  
      `pip install -r requirements.txt`
      
  5. Activating the conda environement.
  
      `conda activate PhyloPy`

# Downloading Executible Dependencies

  1. Click [here](https://drive.google.com/file/d/1T2vvrFE4WY0ViiUXWEtoperxZn3z9DSS/view?usp=sharing) to download the dependencies
  
  2. Extract archive in GitHub Repo dir, should contain execs/ and records.db
  
  3. Change permissions in each executible file to "allow executing as program"

# Cloning the GitHub Repo

  On GitHub, navigate to the main page of the repository.

Above the list of files, click  Code.

"Code" button
To clone the repository using HTTPS, under "Clone with HTTPS", click . To clone the repository using an SSH key, including a certificate issued by your organization's SSH certificate authority, click Use SSH, then click . To clone a repository using GitHub CLI, click Use GitHub CLI, then click .

The clipboard icon for copying the URL to clone a repository
The clipboard icon for copying the URL to clone a repository with GitHub CLI
Open Terminal.

Change the current working directory to the location where you want the cloned directory.

Type git clone, and then paste the URL you copied earlier.

$ git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
Press Enter to create your local clone.

$ git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
> Cloning into `Spoon-Knife`...
> remote: Counting objects: 10, done.
> remote: Compressing objects: 100% (8/8), done.
> remove: Total 10 (delta 1), reused 10 (delta 1)
> Unpacking objects: 100% (10/10), done.
Cloning a repository to GitHub Desktop
On GitHub, navigate to the main page of the repository.

Above the list of files, click  Code.

"Code" button
Click  Open with GitHub Desktop to clone and open the repository with GitHub Desktop.

"Open with GitHub Desktop" button
Follow the prompts in GitHub Desktop to complete the clone.

For more information, see "Cloning a repository from GitHub to GitHub Desktop."

Cloning an empty repository
An empty repository contains no files. It's often made if you don't initialize the repository with a README when creating it.

On GitHub, navigate to the main page of the repository.

To clone your repository using the command line using HTTPS, under "Quick setup", click . To clone the repository using an SSH key, including a certificate issued by your organization's SSH certificate authority, click SSH, then click .

Empty repository clone URL button
Alternatively, to clone your repository in Desktop, click  Set up in Desktop and follow the prompts to complete the clone.

Empty repository clone desktop button
Open Terminal.

Change the current working directory to the location where you want the cloned directory.

Type git clone, and then paste the URL you copied earlier.

$ git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
Press Enter to create your local clone.

$ git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
> Cloning into `Spoon-Knife`...
> remote: Counting objects: 10, done.
> remote: Compressing objects: 100% (8/8), done.
> remove: Total 10 (delta 1), reused 10 (delta 1)
> Unpacking objects: 100% (10/10), done.

# FAQ

  1. Import Error: Cannot import package "SeqMotifFace"
  
   Open terminal:
   
    `conda install -c etetoolkit ete3`
