# File Recovery via RStore Snapshots

## Important Note

- These instructions only apply to RStore (Tier 1) storage.
- Snapshots are taken **daily at 7PM EST**.
- Snapshot retention period is **3 months**. Data deleted after the retention period is **not recoverable**.
- Do not modify or delete any snapshots. There is no snapshot of snapshots.

## Recovery via Windows

1. Open **File Explorer**.
1. Click **This PC**.
1. Click the **...** button and select **Map network drive...**.
1. In the **Folder** text box, enter your RStore path and click **Finish**.
   - Example: `\\rstore.it.tufts.edu\tts_rsch_NAME$`
1. In your mapped drive, create a folder called **RECOVERY**.
1. Click the **...** button and select **Properties**.
1. Click the **Previous Versions** tab.
1. Under **Folder versions**, highlight the date you want to recover data from and click **Open**.
1. Drag and drop the folder or file you want to recover into the **RECOVERY** folder.

## Recovery via Mac OS 14 or 15

1. Open the **Terminal**.

1. Navigate to the `.snapshot` directory for your RStore share:

   ```
   cd /Volumes/tts_rsch_NAME$/.snapshot/
   ```

1. List available snapshots to find the date you need:

   ```
   ls -al
   ```

1. Create a **RECOVERY** folder in your RStore share:

   ```
   mkdir /Volumes/tts_rsch_NAME$/RECOVERY/
   ```

1. Copy your data from the snapshot into the RECOVERY folder.

   Recover a folder recursively:

   ```
   cp -r /Volumes/tts_rsch_NAME$/.snapshot/DATE/FOLDER_NAME /Volumes/tts_rsch_NAME$/RECOVERY/
   ```

   Recover a single file:

   ```
   cp /Volumes/tts_rsch_NAME$/.snapshot/DATE/FILENAME /Volumes/tts_rsch_NAME$/RECOVERY/
   ```
