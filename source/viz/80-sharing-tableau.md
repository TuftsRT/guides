# Sharing and Presenting Tableau Visualizations and Dashboards

Once you have created a table, visualization, or dashboard, what do you do with it? How do you share it with people? This page presents an overview of different methods for using Tableau visualizations in papers, presentations, websites, and more.

## Objectives

- Learn how to bundle your workbooks and all input files into a **Tableau Packaged Workbook** to share your work with others who have access to Tableau Desktop
- Learn how to **export individual visualizations** from worksheets as image files
- Learn how to activate **Presentation Mode** to display Tableau dashboards in a full-screen display suitable for presentations
- Learn how to use Tableau Public for **making dashboards publicly available,** including an option to embed a dashboard in a personal or institutional website

## Using Tableau Packaged Workbooks to Send Tableau Projects to Others with Tableau Desktop

If you want to share your workbook with someone who also has Tableau Desktop installed on their computer, you would also need to send the input data and any other input files used to create the dashboard, such as images. This can be very clumsy for both you and the person you're sending it to.

Luckily, there's an easier way. Tableau has an option to bundle an entire project—including both your Tableau workbook and all the input files—into a single save file, using a **Tableau Packaged Workbook.**

To export your work to a Tableau packaged workbook, click the file menu and select "Export Packaged Workbook..." A window will appear asking you to name your packaged workbook and specify where you want to save it. When ready, click save. You will now have a single file in your save location with the extension .twbx. You can send this file to collaborators without having to send anything else.

For more information on Tableau packaged workbooks, click [here](https://help.tableau.com/current/pro/desktop/en-us/save_savework_packagedworkbooks.htm).

## Exporting Individual Visualizations and Tables

If you would like to use a Tableau visualization in an academic presentation or poster, you'll probably want to convert it into an image file, such as a .png or a .jpeg. There are two ways to do this:

1. Export your worksheet directly to an image file. This is quick and easy, but unfortunately Tableau does not make it easy to control the size and shape of your exported image.
1. Move your worksheet into its own dashboard for more flexible export options.

### Option 1: Exporting directly from a worksheet

To export your worksheet as an image, first navigate to your worksheet by clicking the corresponding tab at the bottom of your screen.

From the menu bar at the top of your screen, click on "Worksheet", then "Export", and then "Image...." A window will pop up asking you which elements of the worksheet you'd like to include in the export, so if you wish to exclude anything, you can simply uncheck it without having to remove that element from the worksheet itself. When you're ready, click save.

Tableau can handle a variety of export formats, including .png, .jpeg, .svg, and .gif. Those familiar with .svg may find it useful when an image needs to be scalable; otherwise, any image format that you are familiar with will work just fine for most cases.

### Option 2: Use a dashboard to control image sizing

For this option, you'll need to create a new dashboard. In the Dashboard panel on the left side of your screen, you sill see a section labeled "Size". Select your desired image dimensions from the drop-down menu, or use the custom height and width controls to specify a size appropriate for your project.

Next, in the Sheets section of the Dashboard panel, choose the worksheet corresponding to the visualization you want to export and drag it to the dashboard. Your visualization will fill the area defined by the dashboard.

Make any additional alterations to your image as needed, and when you're ready, go to the Dashboard menu at the top of your screen. Select "Export Image..." and use your system's default save window to name your exported image, choose your image type, and select a save location.

## Including Dashboards in Presentations

If you're giving a presentation and want to show your dashboard to your audience, Tableau has a convenient full-screen presentation mode. To activate it, locate the following icon in the toolbar at the top of your screen:

<img src="https://tufts.box.com/shared/static/vp1c7iumkbjzno0qwbf3ji3jhdzykqki.png" alt="Presentation mode icon">

Click on this icon to activate presentation mode. You can also press F7 on your keyboard.

To exit presentation mode, press escape or F7. (You could also locate the presentation mode icon, which moves to bottom-right of your screen during presentation mode.)

## Sharing a Tableau Dashboard Publicly on the Internet

Before reading this section, please make careful note of the following:

> **Important!** If you want to make a Tableau dashboard public, there must be no sensitive or protected data not just in your dashboard, but in *any of the input files used to create your dashboard.* This is because publicly available Tableau dashboards can be downloaded *with their source files*, and anyone can see anything that was in the input files used to create the dashboard.

The main vehicle for sharing Tableau dashboards publicly is **Tableau Public:**

https://public.tableau.com/app/discover

Tableau Public includes a server where you can upload your dashboards for public view.

Once your dashboard is uploaded to Tableau Public, Tableau also includes tools for embedding your dashboard on your own personal websites.

### Signing up for Tableau Public

Tableau Public is managed separately from Tableau Desktop licensing and requires a separate account. Tableau Public accounts are also not managed by Tufts, so you will need to create an individual account directly through the Tableau Public website. This may be your own personal account or an account to be shared with a laboratory or research group.

To create an account, navigate to the URL linked above, and click on "Sign Up for Tableau Public". Feel free to use your Tufts email if you only plan to use this account for Tufts-related work. Otherwise, you may choose to use your personal email address if you prefer.

### Uploading your Dashboard to Tableau Public

> Note: for convenience, it is recommended that you export to a Tableau packaged workbook before this step. Follow the steps in the "Using Tableau Packaged Workbooks..." section above, then close Tableau and open the Tableau packaged workbook version of your project that you just saved.

To upload your project to Tableau Public, you will first need to open your project in Tableau Desktop. Navigate to the Data Source tab, and look to the top right-hand corner of the screen for the word "Connection", below which there are two check boxes that say "Live" and "Extract". By default, "Live" is checked, but you need to change this to "Extract" in order for Tableau to create an "extracted" copy of your data that it can carry with it to Tableau Public.

Next, navigate to the dashboard you want to upload. From the "Server" menu, select "Tableau Public" and then choose "Save to Tableau Public" from the submenu. You will then be prompted to sign in to Tableau Public. Be sure to use the separate log-in credentials you created for Tableau Public (See "Signing up for Tableau Public" above).

Depending on your version of Tableau Desktop, you may or may not get sent directly to the Tableau Public website. If you're using an older version of Tableau, you may get a message that says: "The workbook was published successfully. Upgrade to the latest version of Tableau Desktop to proceed directly to a viz after publishing." If this is the case, simply navigate to the Tableau Public website in your web browser at [public.tableau.com](public.tableau.com) and sign in using your Tableau Public credentials.

To view your dashboard on Tableau public, click the account icon on the top right-hand corner of the website:

<img src="https://tufts.box.com/shared/static/zxzf7c8adxrlck1oq5sw2813ugdaxn1q.png" alt="Account Icon in Tableau Public">

A drop-down menu will appear. Click "My Profile" and you will be directed to your personal page on Tableau Public. You should see the dashboard you created under "Vizzes".

### Sharing a Dashboard that has been Uploaded to Tableau Public

From your profile in Tableau Public, you should now click on the Dashboard that you want to share.

In the dashboard view, you will see a row of icons at the top-right corner. Look for the Share icon, which looks like this:

<img src=https://tufts.box.com/shared/static/ilks2ewi5s81fen65frl53mk6ahsmlci.png, alt= "Share Icon in Tableau Public">

Click on the share icon and a pop-up window will appear with two share options:

- A shareable web link with a "Copy Link" button next to it
- An "Embed Code" option, with a "Copy Embed Code"
  - Use this option if you have have a website where you would like to embed your dashboard. You can copy and paste the html code directly into the body of your website in an html or markup editor where you edit your website. We have done this below; note how this creates a fully interactive copy of your dashboard in the body of your website.

<div>
<div class='tableauPlaceholder' id='viz1747940314047' style='position: relative'><noscript><a href='#'><img alt='Customer Satisfaction Dashboard ' src='https:&#47;&#47;public.tableau.com&#47;static&#47;images&#47;Tu&#47;TuftsTableauTutorial-JumbosCroissants&#47;CustomerSatisfactionDashboard&#47;1_rss.png' style='border: none' /></a></noscript><object class='tableauViz'  style='display:none;'><param name='host_url' value='https%3A%2F%2Fpublic.tableau.com%2F' /> <param name='embed_code_version' value='3' /> <param name='site_root' value='' /><param name='name' value='TuftsTableauTutorial-JumbosCroissants&#47;CustomerSatisfactionDashboard' /><param name='tabs' value='no' /><param name='toolbar' value='yes' /><param name='static_image' value='https:&#47;&#47;public.tableau.com&#47;static&#47;images&#47;Tu&#47;TuftsTableauTutorial-JumbosCroissants&#47;CustomerSatisfactionDashboard&#47;1.png' /> <param name='animate_transition' value='yes' /><param name='display_static_image' value='yes' /><param name='display_spinner' value='yes' /><param name='display_overlay' value='yes' /><param name='display_count' value='yes' /><param name='language' value='en-US' /></object></div>                <script type='text/javascript'>                    var divElement = document.getElementById('viz1747940314047');                    var vizElement = divElement.getElementsByTagName('object')[0];                    if ( divElement.offsetWidth > 800 ) { vizElement.style.width='1000px';vizElement.style.height='827px';} else if ( divElement.offsetWidth > 500 ) { vizElement.style.width='1000px';vizElement.style.height='827px';} else { vizElement.style.width='100%';vizElement.style.height='927px';}                     var scriptElement = document.createElement('script');                    scriptElement.src = 'https://public.tableau.com/javascripts/api/viz_v1.js';                    vizElement.parentNode.insertBefore(scriptElement, vizElement);                </script>
</div>
