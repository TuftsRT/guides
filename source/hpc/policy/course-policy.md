# HPC Course Policy

## Purpose and Scope

The HPC infrastructure is a critical tool in enhancing educational experiences, providing students and faculty access to advanced computing resources for research, learning, and academic development. The purpose of this policy is to outline the guidelines and standard procedures for the use of High-Performance Computing (HPC) resources and services in support of teaching and learning activities at Tufts University.

This policy applies to all users of HPC resources, including students, faculty, staff, and other authorized personnel.

## Notable changes

For previous users of the Tufts cluster these are the items that deviate from previous practice the most. They are all covered below, but we want to highlight these differences.

- Courses will be given a unique, new project folder each semester. At the end of the semester this folder will be archived.
- Timeline, we are requesting all instructors request access for their course 4 weeks prior to the start of the semester.
- Course project folders from before 7/1/2025 will be archived on Tier 2 for 12 months. RT will contact each instructor with additional details.

## Guidelines

This policy starts 7/1/2025, however we will be flexible with existing instructor needs as we implement this. Please contact RT with any concerns.

**Each semester, instructors are required to submit a new request to use HPC resources for their course using the
[HPC Course Request Form](https://tufts.qualtrics.com/jfe/form/SV_d7o0UZFgK1PFXnv):**

### Timeline

- HPC course setup request needs to received by RT 4 weeks prior to the first day of class.
- Full student list must be provided at least 1 week prior to the first day of class.
- Course project folders will be archived 30 days after the end of the semester.

### Capacity

- RT can currently support 10 non-concurrent courses using the cluster
- Research technology will try our best to accommodate course HPC resource needs. When in high demand, HPC course requests are prioritized based on order of request.

#### Compute

Each course may request a pool of computing resources be reserved for in-class/lab sessions during the semester.\
Public partition HPC computing resources cannot be reserved for use outside of the scheduled class time, or for
courses that do not meet synchronously.

Maximum reserved resources per course:

- 128 CPU cores
- 1TB CPU memory
- 20 GPUs (16GB GPU Device VRAM each)
- For any additional computational resource needs:
  - CPU resources - please reach out to RT for further information on options.
  - GPU resources - hardware purchases supported by school, department, or instructor may be needed.

#### Storage

Each course may request up to 1TB of Tier 1 HPC storage, which will be setup as follows:

- A clean project folder will be created each semester for each course
  - Within the course folder, a "shared" subfolder will be created for any files need to be accessed by the entire class (read-only to students)
  - Within the course folder, each student will have their own subfolder which is only accessible by this student, the instructor, and course TA(s).
  - Within the course folder, group project folders can be created upon additional request from the instructor.
  - This course storage quota does not count towards the instructors research storage quota for the charge back model.
- The course project folder will be archived 30 days after the end of semester
  - Moved to Tier 2 Research - retain one year, then deleted

### Software

- The RT team provides its best effort to have requested compatible open-source software or packages installed and set up 2 weeks before first day of class (under the condition the request was received at least **4 weeks** prior to first day of class.)
- Any changes to a course software environment setup must be communicated to RT 2 weeks prior to the day of in-class
  session.
- Instructor and/or designated TA(s) will be responsible for testing software setup prior to the day of in-class session.
- New commercial software needs to be purchased by instructor or department at least **4 months** before first day of class.
- The general practice is to not change or update the versions of software in use for courses during the semester

### Support

Support for students using the HPC cluster for coursework is provided via a shared model between the course Instructor, TA's and RT staff.

- Instructors or TA's provide assistance with:

  - Instructing students on the use of the software
  - Debugging student code or assignment problems
  - Addressing any problems that arise during class period, including initial debugging
  - Acting as a single point of contact with RT about any common issues students encounter on the HPC cluster, and
    communicate the solutions RT provides to students.

- RT provides assistance with:

  - Setting up student accounts
  - Course HPC environment setup (storage, software installation)
  - The in-class "Intro to HPC" sessions
  - Troubleshooting HPC system wide problems as reported by the Instructor, typically indicated by multiple students
    having the same problem

- For Capstone or Project-Only courses, all support requests should be made by a designated student or team leader
  who serves as a single point of contact.

- In-class Sessions

  - A Virtual "Intro to HPC" in-class session is available. Requests need to be received by RT 4 weeks prior to the scheduled in-class session.
    - When scheduling permits, in-person sessions may also be available.
  - Please contact tts-research@tufts.edu about special topic sessions.

### Availability

RT makes every effort to provide a reliable platform for research and teaching on the Tufts HPC Cluster. However the nature of HPC means that the system experiences failures, or degraded performance at times. End user support and system administration is provided during normal business hours.

- Courses running on the HPC cluster are not immune to scheduled or unscheduled cluster-wide outages.
- Courses are effected by the regular HPC maintenance period, and annual MGHPCC shutdown.
