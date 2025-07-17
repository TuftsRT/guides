# Introduction to Dashboards

You've now created some bar charts and summary tables, but if you're looking to put them on a website, you'll probably want to bring them all together in a single, attractive layout. You may even want to include interactive elements that allow your user to select what data they want to look at. 

For this, we'll need to create a **dashboard.** 

## Objective

- Arrange the contents of multiple worksheets into a visually pleasing dashboard 
- Include an interactive filter for your end-user to select subsets of the data they're most interested in seeing


## Creating a New Dashboard

To create a new dashboard, locate the "new dashboard" button at the bottom of your screen, right next to the "new worksheet" button you used in the previous lessons. It will look like this:

<img src="https://tufts.box.com/shared/static/wr4fne9nzgn7mswdi64ribndioenx1ad.png" alt="New Dashboard button">

You can check that you've found the correct button by hovering your mouse over it. If you're in the right place, some text will pop up that says "New Dashboard". 

Once you open a new Dashboard, your screen will look like this:

<img src="https://tufts.box.com/shared/static/izmj9krfw9fienpgwqdc93ykswdupg0e.png" alt="New Dashboard Screen">

## Adding Worksheets to our Dashboard

On the left of your screen, you'll see the **Dashboard Pane**, with subheadings for Size, Sheets, and Objects. We'll stick with the default dashboard size for now, which is optimized for dashboards intended to be viewed on a Desktop Browser. 

Now take a look the Sheets section, where you'll see a list of all the worksheets we've created so far. These are the building blocks of our dashboard; dashboards are a place to assemble, arrange, and present the information you created in worksheets.

Let's start with the "Counts by Gender" table we created in the lesson on tables.

### Adding our Counts by Gender Table

Adding a the contents of a worksheet to a dashboard is a breeze. Simply click on "Counts by Gender" in the Sheets section of the Dashboard pane, and drag it to the view panel, i.e. the place where it says "Drop sheets here". 

You should now see the table appear in your dashboard.

### Adding our Counts by Rating Bar Chart

A useful feature of dashboards is that they can be interactive, allowing the user to choose the data they want want to focus on. 

It would be nice to include some interactive features in our dashboard too. For example, it would be nice if we could give the user a way to filter the responses to a specific gender, so they can see if all demographics are being served equally. Let's make a quick detour back to our "counts by rating" chart and add one additional element to ensure this kind of interactivity will be enabled when we add it to our dashboard.

#### Adding an interactive filter

Navigate back to the "Counts by Rating Chart" by clicking on the corresponding worksheet tab at the bottom of your screen. Once there, find "Gender" in the data pane, and drag it to the Filters shelf. 

A pop-up window will appear:

<img src="https://tufts.box.com/shared/static/8w5rsgsho226h1z89y4lntomoya7p207.png" alt="Gender filter pop-up window">

Make sure all options (Male, Female, and Non-Binary) are checked, and click "OK". Your chart will still look the same as before, but you'll notice that there's now a pill for the gender field in the Filters shelf. 

We now have one more step to ensure that this filter appears on our dashboard as an interactive element. Right-click the Gender pill in the filters shelf, and from the pop-up menu, select "Show Filter", as shown below:

<img src="https://tufts.box.com/shared/static/lcadwyzfsnbwp5t1rgad5jjcpuo82sb8.png" alt="Selecting the Show Filter option">

Once you click that, you should notice a new panel of check boxex appear next to your chart in the worksheet view, like this: 

<img src="https://tufts.box.com/shared/static/gbzfy82b0l6wt3xbzkui46ibwupm76tc.png" alt="Response counts by rating chart with gender check boxes">

These check boxes will now follow our chart into our dashboard and will be available to our end-user to interact with.

Now we're ready to add our interactive chart to our dashboard!

#### Placing the Response Count Bar Chart on the Dashboad

Navigate back to your dashboard by selecting the corresponding tab at the bottom of your screen, and locate the Counts by Rating Chart in the sheets section of your Dashboard Pane.

In a moment, we'll drag and drop this chart into our dashboard, just as we did with the Counts by Gender table. This time, however, where we place the chart is going to matter. Since we already have one element in our dashboard (the Counts by Gender table that we added previously), we need to tell Tableau were to put the new chart in relation to this existing element.

Click and drag the Counts by Rating chart to the dashboard, but don't release your finger from the mouse just yet. Watch what happens when you move your cursor to different parts of the screen. You'll see different shaded regions appear, corresponding to different placement options for your chart. 

For our dashboard, we want our Counts by Rating chart to be on the left side of the dashboard, so you should release your finger from the mouse button when your screen looks like this:

<img src="https://tufts.box.com/shared/static/igkuv2kz9hng7cx9ho7dwgxr02kjbaez.png">

If you did it correctly, your dashboard should now look like this:

<img src="https://tufts.box.com/shared/static/78kkwinizjvrvj49wxozl9dbaz3vxym3.png" alt="Dashboard after placing the Counts by Rating chart">

#### Rearranging the Dashboard Elements

Let's try to get this looking better by rearranging our chart, table, and interactive filter in a more visually appealing way. We can start by moving the filter closer to the chart.

Click on the headading of the interactive filter to select the entire filter element. When selected, it should look like this:

<img src="https://tufts.box.com/shared/static/do3sgqod7a1ib7oa5ix23ow9nwvgwisc.png" alt="Interactive filter selected">

When selected, you can move it by clicking and dragging here: 

<img src="https://tufts.box.com/shared/static/sfn8b7fxm6q2li3grbzlkgu20sawhdum.png" alt="Click and drag the filter here">

Guides will appear to help you place it correctly. Follow the example below to move it below your Counts by Gender table, next to your Number of Responses by Rating bar chart. 

<img src="https://tufts.box.com/shared/static/4xqgh96n4yzr9zpesenr6fn6dldvx5ya.png" alt="Placing the filter next to the chart and below the gender counts table">

Your dashboard will now look like this:

<img src="https://tufts.box.com/shared/static/qiifb5vupdpfw659wz3c1o7hjao2k3h3.png" alt="Dashboard after moving the filter">

Now let's resize a few elements to make this look better:
- Click and drag the tob border of the filter to resize it so that you can see the full filter. 
- Select the bar chart, and widen it by clicking and dragging the right border of the chart. 

With the resized elements, your dashboard will be looking more and more like a professional layout!

<img src="https://tufts.box.com/shared/static/butphvad08qp2oeqlqgmes40s6wnbf8o.png" alt="Dashboard with resized elements">

### Testing out the Interactivity

Now take a moment to test out the interactivity of our dashboard. Try checking and unchecking some of the boxes in the Gender filter. Notice how the chart adjusts to show you data for only the genders you've selected? 

If you didn't have this interactive filter enabled, you'd have to create a separate chart for every gender, which would add a lot of unnattractive complexity to the dashboard. Interactive elements help simplify your visual layout while still presenting a large amount of information. They can also engage your end-user by making them an active partner in exploring your data. 

## Adding a Title to the Dashboard

Now let's add a text element to our dashboard that can serve as a title or header. 

In the Dashboard pane at the left side of your screen, locate the section that says "Objects" and find within it the text object. 

<img src="https://tufts.box.com/shared/static/10jih4frav0drc6euivfaqkb813z31km.png" alt="The Text object within the Objects pane">

Drag it to your dashboard, using the guides to place it at the top of your screen, like so:

<img src="https://tufts.box.com/shared/static/igan7odcjr7xyhcjiebilqxmd8ujk9ux.png" alt="Drag the text object to the top of the dashboard">

A pop-up window will appear asking you to specify the text for your text object. In this window, type the following:
- In size 48 font, type "Jumbo's Croissants"
- In size 16 font, on a new line, type "Customer Satisfaction on a Scale from 1 (Poor) to 10 (Outstanding)"

<img src="https://tufts.box.com/shared/static/qw1mwrl4m54fk3g8ssoe1nw7baopmgge.png" alt="Edit title text">

Then press OK to see your text appear in the text object on your dashboard. 

You'll probably notice that the text object you just created is a bit too big. That's easy to fix! Select the text object, just as before, and click and drag the bottom border to resize it. 

Your dashboard now has a title:

<img src="https://tufts.box.com/shared/static/uwo43lt08irxls3h4c9mbsgv8euyzc5y.png" alt="Dashboard with text object added for title">

## Adding an Image

Jumbo's Croissants has a logo, and Jumbo would love to show it off. 

<img src="https://tufts.box.com/shared/static/ijl45unwq3m1y0nl1nhn7gfe4bbjrppm.png" alt="Jumbo's Croissants Logo">

You can access this image at [this link.](https://tufts.box.com/shared/static/ijl45unwq3m1y0nl1nhn7gfe4bbjrppm.png) Once you click on the link, right-click on the image and select "Save image as" to save it to your computer. 

Adding an image object works exactly the same way as adding a text object. In the objects section of the Dashboard pane, click on "Image" and drag it to the upper left-hand corner of your dashboard, as shown here:

<img src="https://tufts.box.com/shared/static/vleyki9ajlpspzzj8yo18t02esfpuhl5.png" alt="Drag image to upper left corner">

A pop-up window will appear from which you can select your image file. Find where it says "Choose and image file..." and click choose, and then navigate to the location on your computer where you saved your image file. When you're ready, click "OK", and the image will appear in the image box you just created.

<img src="https://tufts.box.com/shared/static/o2u7zloxi1fraar9x74wmw4hxgeup8sh.png" alt="Dashboard with image added, before the image is resized">

When the image first loads, it appears on your dashboard the same size as it is in the image file. To make it scale to the size our image object, right-click on the image and select "Fit Image" from the pop-up menu.  

<img src="https://tufts.box.com/shared/static/nilr5n8nn04qf3pl3bvzlzryzdus54pz.png" alt="Right-click on the image and select Fit Image">

You will now see the full image. Let's also resize the image object by clicking and dragging the right border and making it smaller. 

<img src="https://tufts.box.com/shared/static/nilr5n8nn04qf3pl3bvzlzryzdus54pz.png" alt="Dashboard with resized and fitted image">

## Final Touches

We're almost done! Let's just make a few additional changes to make things look more professional and visually appealing.

### Reducing Visual Clutter

Removing unnecessary elements can make for a more streamlined and professional-looking dashboard. For example, when we brought in our bar chart, it came with the title that we gave it in our worksheet. This is the text that says "Number of Responses by Rating" at the top of our chart. This is probably unnecessary for our end-user to understand the chart, since our Dashboard already tells them that this is a customer satisfaction survey on a 1-10 scale, and the y-axis already indicates that we're charting the number of responses. So let's get rid of the chart title to reduce visual clutter. 

To get rid of the chart title, simply right-click on where it says "Number of Responses by Rating" and select "Hide Title". 

<img src="https://tufts.box.com/shared/static/jev3mlpj888bqn4qq1t7kr287mpwhzcl.png" alt="Hide the bar chart title">

Another place where we can get rid of unnecessary text is in the counts by gender table. Notice the word "Gender" is still hanging out at the top of the column that lists the different genders. But given that we already have a descriptive title for our chart (i.e. "Charts by Gender"), and the fact that it's already clear what the table is showing, we really don't need this extra instance of the word "Gender". 

To remove it, right-click it and select "Hide field labels for rows", as shown in the image:

<img src="https://tufts.box.com/shared/static/0wfxgolv7q39jqg7qsp97svfe12nrf0y.png">

### Improving the Filter Title

We can also make a few improvements to the title of the interactive filter in the bottom right of your dashboard. At the moment, it just says "Gender", but let's change it to "Select Gender" to draw the end-user's attention to it as an interactive element. 

Right-click the title (where it says "Gender", just above the check boxes) and select "Edit Title" to get a pop-up window where you can edit the text for the interactive filter title. 

<img src="https://tufts.box.com/shared/static/ify0xd5wng87gw23hnp3hj3mraw4pmom.png" >

Change the text to "Select Gender", and then change the font size to 16. You can also change the font to "Tableau Light" to be visually consistent with the rest of the dashboard text. When you're ready, click "OK". 

### Rename the Dashboard

As a final step, let's give our dashboard a name. At the bottom of the screen, locate the tab for your dashboard. It should currently say "Dashboard 1". Just as you did with worksheets, you should right-click on the tab and select rename. Rename it "Customer Satisfaction Dashboard", or whatever name feels best to you. 


## Admire your Completed Dashboard

And we're done! It's time to admire your completed dashboard. 

<img src="https://tufts.box.com/shared/static/ueisbazl30vp719852l8vacb4sx1c2hi.png" alt="Completed Dashboard">


## Next Steps

We've created some tables and visualizations, but how do we actually use them and share them elsewhere? The next lesson will cover sharing and presenting your creations in Tableau, as well as embedding Tableau Dashboards in websites.

## Additional Resources

- [Tableau Desktop and Web Authoring Help (Official Tableau documentation))](https://help.tableau.com/current/pro/desktop/en-us/default.htm)