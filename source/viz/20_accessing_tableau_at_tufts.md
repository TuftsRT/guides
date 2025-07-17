# Accessing Tableau at Tufts

## Objective
Learn how to access Tableau Desktop at Tufts, whether you are an administrator, student, faculty member, or research staff.



## Understanding Tableau Public vs. Tableau Desktop vs. Tableau Desktop Public Edition
- **Tableau Public**: A free web-based version of Tableau where all published data is publicly accessible. This is not recommended for most university work as it does not ensure data confidentiality. However, if you would like to make your data publicly available and it is suitable for public access, you can publish your dashboard to the Tableau Public Gallery.
- **Tableau Desktop**: A downloadable desktop version of Tableau that allows private data analysis and visualization. This is the recommended tool for all university-related work, including research, teaching, and internal business purposes. It can be accessed in-person via the Tufts computer labs or remotely through the Tufts remote computing insfrastructure (see below for more information).
- **Tableau Desktop Public Edition**: This is a slightly more limited version of Tableau Desktop that's free for students and researchers. You can download this and use it with an academic license, in which case it you may only use it for academic purposes. It allows you to work privately and upload your dashboards to Tableau Public, while limiting the options for sharing that are available in the full Desktop edition. There are some other minor differences as well; [see here for a full comparison between Tableau Desktop and Tableau Desktop Public Edition](https://help.tableau.com/current/pro/desktop/en-us/desktop_comparison.htm?_gl=1*jmvnxr*_ga*MTI5NTYxOTU4Ni4xNzQzMTY5Mjcw*_ga_8YLN0SNXVS*czE3NTI2NzQ4NDAkbzEzJGcxJHQxNzUyNjc0ODcxJGoyOSRsMCRoMA..*_gcl_au*MTQzMjYwNTk2MS4xNzUxODk1Mzgz). This is a good option if you want to learn Tableau and don't need to share your workbooks with others, or if you plan to share your workbooks publicly via Tableau Public. Note that you can still show your dashboards during presentations and team meetings using your own Tableau Desktop Public Edition license.

The walkthrough below will use Tableau Desktop. You can see instructions on how to use Tableau Public here. 



## Ways to Access Tableau Desktop
Your method of accessing Tableau Desktop depends on your role at Tufts and your intended use:



### Administrators or Staff Using Tableau for Internal Business Purposes

Using Tableau for internal business purposes requires a paid Tableau license. This is managed separately from Tableau used for academic purposes (research or coursework). Administrators should not follow the instructions for academic users below. Instead, please go to: https://access.tufts.edu/tableau


### Students, Faculty, and Research Staff Using Tableau for Academic Research or Coursework
Students, faculty, and research staff have three primary options for accessing Tableau Desktop:

#### **Option 1: Download and Install Tableau Desktop Public Edition**
- **Overview** Tableau Desktop Public Edition may be suitable for you if you don't mind the restricted sharing options and plan to use it solely for academic and research purposes. You can use it to practice and learn Tableau on your local computer, create dashboards that are suitable for public sharing via Tableau Public, or present your findings in-person to your collaborators while logged in to your own personal account. However, note that Tableau prohibits using your student license for work related to outside jobs or internships.
- **Eligibility**: Free for students and educators with a valid `.edu` email address. 
- **Steps**:
  1. Visit the [Tableau Academic Program website](https://www.tableau.com/community/academic).
  2. Request a free license and download Tableau Desktop Public Edition.
  3. Install Tableau Desktop Public Edition on your computer and activate it with the license key provided.


#### **Option 2: Access Through Tufts Remote Computing Infrastructure**
- **Overview**: This option allows you to use a Tufts-managed Windows computer remotely from your personal computer, which in turn gives you full access to Tableau Desktop.
- **Steps**:
  1. Install the Omnissa (formerly VMware) Horizon desktop client (recommended) or use the web application. For instructions, visit the [Tufts Remote Computing page](https://access.tufts.edu/remote-computing-infrastructure).
  2. Log in to the Tufts Remote Computing Infrastructure using your Tufts credentials.
  3. Launch Tableau Desktop from the remote environment. Tableau can be found by searching for "Tableau" in the Windows search bar or by clicking on the Tableau desktop shortcut.
        - **NOTE:** Virtual sessions may take some time to fully load. If Tableau doesn’t immediately appear in the Windows search bar or on the Desktop, wait and try again. It may take up to a minute for the system to index files and make Tableau visible.
 - **Accessing and Storing Files**: 
    - Do not attempt to store files or save your work in the remote computer's file system. Files will be deleted after the session ends. Instead, consider either of the following two options for easy and convenient access to your files:
        - **Cloud Storage with Box**
            - Box is available to all Tufts Faculty, Staff, and Students.
            - When you install Box on your machine, it functions like any other folder in your file system. In most cases, you can securely save your data and your work here. However, if you are working with sensitive or individually identifiable data there may be other considerations. For more information, see the Tufts [Data Storage Finder.](https://access.tufts.edu/data-finder)
            - Files can then be accessed from any other computer by logging into Box either through a web browser or a desktop client. The Box desktop client is already installed on the remote computer.
            - For more information and installation instructions, click [here](https://access.tufts.edu/box). 
        - **Allow the Remote Computer to Access Local Files**
            - For those using the VMware Horizon desktop client, you can configure the settings to grant the remote computing environment to have direct access to files saved on your local computer. This requires a one-time change to the settings and then you will not have to repeat this step.
            - To enable this:
                1. Open VMware Horizon Client and log in to vdi.it.tufts.edu.  
                2. After logging in, open "Settings" from the navigation panel on the left-hand side of the window. 
                3. In Settings, click on "Applications". Locate the setting "open local files in hosted applications" and set it to "on".
                4. [Optional] If you have a multi-monitor setup, you may wish to restrict the remote session to a single monitor so that you can still interact with your own computer on the other monitors. To do this, open settings and from the menu on the left select the remote resource you're using (e.g. "Engineering Lab"). Click on it and you should see a drop-down menu labeled "Display". Select "Fullscreen - Single Monitor" from this menu. 
            - Once you have enabled the remote computer's access to your local file system, log in to the remote computer and open file explorer. Click on "This PC" and select "Network Drive (Z:)". This will take you to your local computer's file system.
            - NOTE: Depending on the size and complexity of your project, as well as the strength of your network connection, working out of a mounted drive may cause your Tableau session to run slowly. If this is the case, you might consider moving your project to a temporary folder on the remote desktop. This can sometimes lead to big improvements in the speed of your Tableau session, but be careful to create frequent back-up copies in your local computer's files, because if your session gets disconnected, you will lose any files that are saved on the remote computer.



#### **Option 3: In-Person Use of the Tufts Data Lab and Other Computer Labs**
- **Overview**: Tableau Desktop is pre-installed on computers at various labs including the [Medford and Boston Tufts Data Labs](https://sites.tufts.edu/datalab/). 
- **Best Use**: Ideal for those who prefer working on-site with ready-to-use Tableau installations.




## Next Steps
With Tableau installed and ready to use, we can now proceed to the next lesson on loading and viewing your data.

## Additional Resources
- [Tableau Academic Program](https://www.tableau.com/community/academic)
- [Tufts Remote Computing Infrastructure](https://access.tufts.edu/remote-computing-infrastructure)
