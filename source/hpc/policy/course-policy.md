# HPC Course Policy
### Purpose and Scope

The HPC infrastructure is a critical tool in enhancing educational experiences, providing students and faculty access to advanced computing resources for research, learning, and academic development. The purpose of this policy is to outline the guidelines and standard procedures for the use of High-Performance Computing (HPC) resources and services in support of teaching and learning activities at Tufts University.

This policy applies to all users of HPC resources, including students, faculty, staff, and other authorized personnel.

### Objectives

To better support teaching and learning and improve resource utilization on Tufts HPC cluster, TTS Research Technology (RT) is announcing the following guidelines to better assist both instructors and RT team with HPC course planning, setup, and support. It will also help RT team manage HPC cluster computational and storage resources and provide prompt support more effectively.

### Noteable changes

For previous users of the Tufts cluster these are the items that deviate from previous practive the most.  They are all covered below, but we want to highlight these differences.

* Courses will be given a unique, new project folder each semester.  At the end of the semester this folder will be archived.
* Timeline, we are requesting all instructors request access for their course 4 weeks prior to the start of the semester.

### Guidelines

Starting 7/1/2025:

**Existing course storage and setup:**

- RT will clean up existing course storage and setup:
  - Archive course folder to Tier 2 storage (inaccessible from cluster), and delete after 12 months.
  - Delete all student folders and data from Tier 1 course storage.
  - Remove all student access to course storage.

**Each semester, instructors are required to submit new requests to use HPC resources for their course :**

- HPC course setup request needs to received by RT 4 weeks prior to the first day of class.

  - Full student list must be provided at least 1 week prior to the first day of class.
  - Research technology will try our best to accommodate course HPC resource needs. When in high demand, all HPC course resources are on "first come first serve" basis.
  - Up to 10 non-concurrent courses
  - [Link to HPC Course Request](#)

- Each course may request a pool of computing resources be reserved for in-class/lab sessions during the semester.  Public partition HPC computing resources cannot be reserved for use outside of the scheduled class time, or for courses that do not synchronously.

  Per course reservation maximums:
  - 128 CPU cores
  - 1TB CPU memory
  - 20 GPUs (16GB GPU Device VRAM each)
  - For any additional computational resource needs:
    - CPU resources - please reach out to RT for further information on options.
    - GPU resources - hardware purchases supported by school, department, or instructor may be needed.

- Each course may request up to 1TB of Tier 1 HPC storage, which will be:

  - A clean project folder will be created each semester for each course
    - Within the course folder, a "shared" subfolder will be created for any files need to be accessed by the entire class (read-only to students)
    - Within the course folder, each student will have their own subfolder which is only accessible by this student, the instructor, and course TA(s).
    - Within the course folder, group project folders can be created upon additinal request from the instructor.
    - This course storage quota does not count towards the instructors research storage quota for the charge back model.
  - The course project folder will be archived 30 days after the end of semester
    - Move to Tier 2 - retain one year, then delete - RT

- Software

  - The RT team provides its best effort to have requested compatible open-source software or packages installed and set up 2 weeks before first day of class (under the condition the request was received at least **4 weeks** prior to first day of class.)
  - Any changes to course software environment setup must be communicated to RT 2 weeks prior to the day of in-class session.
  - Instructor and/or designated TA(s) will be responsible for testing software setup prior to the day of in-class session.
  - New commercial software needs to be purchased by instructor or department at least **4 months** before first day of class.
  - The general practice is to not change or update the versions of software in use for courses during the semester

- In-class Sessions

  - A Virtual "Intro to HPC" in-class session is available. Requests need to be received by RT 4 weeks prior to the scheduled in-class session.
    - When scheduling permits, in-person sessions may also be available.
  - Please contact tts-research@tufts.edu about special topic sessions.

- Support
    Support for students using the HPC cluster for coursework is provided via a shared model between the course Instructor, TA's and RT staff.  

  - RT provides assistance with: access and storage setup, class environment setup, in-class intro HPC sessions (remote), system debugging, .etc.
  - During course period, instructor and course TAs are responsible to collect and attempt the initial debugging when students' issues arise.
  - Instructors and/or designated TAs should be in contact with RT about large scaled common issues students run into on HPC cluster, and communicate the solutions RT provides to students.
  - For Capstone or Project-Only courses, all support requests should be made by a single designated student or team leader.

- Availability
    RT makes every effort to provide a reliable platform for research and teaching on the Tufts HPC Cluster.  However the nature of HPC means that the system experiences failures, or degraded performance at times.  End user support and system administration if provided during normal business hours.
  - Courses running on the HPC cluster are not immune to scheduled or unscheduled cluster-wide outages.
  - Courses are effected by the regular HPC maintenance period, and annual MGHPCC shutdown. 
