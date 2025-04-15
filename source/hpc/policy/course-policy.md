# HPC Course Policy
### Purpose and Scope

The HPC infrastructure is a critical tool in enhancing educational experiences, providing students and faculty access to advanced computing resources for research, learning, and academic development. The purpose of this policy is to outline the standardized guidelines and procedures for the use of High-Performance Computing (HPC) resources and services in support of teaching and learning activities at Tufts University.

This policy applies to all users of HPC resources, including students, faculty, staff, and other authorized personnel.

### Objectives

To better support teaching and learning and improve resource utilization on Tufts HPC cluster, TTS Research Technology (RT) is announcing the following guidelines to better assist both instructors and RT team with HPC course planning, setup, and support. It will also help RT team manage HPC cluster computational and storage resources and provide prompt support more effectively.

### Guidelines

Starting MM/DD/YYYY:

**Existing course storage and setup:**

- RT will clean up existing course storage and setup:
  - Archive course folder to Tier 2 storage (inaccessible from cluster), and delete after 12 months.
  - Delete all student folders and data from Tier 1 course storage.
  - Remove all student access to course storage.

**Each semester, instructors are required to submit new requests for HPC resources and setup:**

- HPC course setup request needs to received by RT 4 weeks prior to the first day of class.

  - Full student list must be provided at least 1 week prior to the first day of class.
  - Research technology will try our best to accommodate course HPC resource needs. When in high demand, all HPC course resources are on "first come first serve" basis.
  - Up to 10 non-concurrent courses
  - [Link to HPC Course Request](#)

- Each course may request a set amount of computing resources for in-class/lab sessions only through out the semester.

  Up to a total of:

  - 128 CPU cores
  - 1TB CPU memory
  - 20 GPUs (16GB GPU Device VRAM each)
  - For any additional computational resource demand:
    - CPU resources - please reach out to RT for further information.
    - GPU resources - hardware purchases supported by school, department, or instructor may be needed.

- No HPC computing resources can be reserved for project-only courses.

- Each course may request up to 1TB HPC Tier 1 (Hot) storage, which will be:

  - Created from scratch for each semester upon request
    - Within the course folder, a "shared" subfolder will be created for any files need to be accessed by the entire class (read-only to students)
    - Within the course folder, each student will have their own subfolder which is only accessible by this student, the instructor, and course TA(s).
    - Within the course folder, group project folders can be created upon additinal request from the instructor.
    - Course storage quota does not count against instructor research storage quota for the charge back model.
  - Cleaned up 30 days after end of semester
    - Move to Tier 2 - retain one year, then delete - RT

- Software

  - RT team will offer best effort to have requested compatible open-source software or packages installed and set up 2 weeks before first day of class (under the condition the request was received at least **4 weeks** prior to first day of class.)
  - Any changes to course software environment setup must be communicated to RT 2 weeks prior to the day of in-class session.
  - Instructor and/or designated TA(s) will be responsible for testing software setup prior to the day of in-class session.
  - New commercial software needs to be purchased by instructor or department at least **4 months** before first day of class.

- In-class Sessions

  - Virtual "Intro to HPC" in-class session is available and request needs to be received by RT 4 weeks prior to the scheduled in-class session.
    - When schedule permits, in-person sessions can be scheduled in advance.
  - Please contact tts-research@tufts.edu about special toptic sessions.

- Support

  - RT offers assistance in: access and storage setup, class environment setup, in-class intro HPC sessions (remote), system debugging, .etc.
  - During course period, instructor and course TAs are responsible to collect and attempt the initial debugging when students' issues arise.
  - Instructors and/or designated TAs should be in contact with RT about large scaled common issues students run into on HPC cluster, and communicate the solutions RT provides to students.
  - HPC course setup is not immune to cluster-wide outage (scheduled or unscheduled).
  - For project-only courses, designated student or team leader will be responsible for the communication with RT team regarding support issues.
