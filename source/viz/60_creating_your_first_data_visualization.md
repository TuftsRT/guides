# Creating Your First Data Visualizations

We're ready to make our first visualizations! In this lesson we'll learn the basics of creating bar charts using both measures and dimensions.

## Objective

- Learn the basics of creating data visualizations in Tableau by creating simple bar charts


## Open a New Worksheet

To create a new data visualization, we first need to open a new worksheet.

In the tray at the bottom left of your screen, locate this new worksheet button:

<img src="https://tufts.box.com/shared/static/fiz7k0izjgomih92uodqzjrs2ss32ovn.png" alt = "Create new worksheet button">

Click this button once to open a new worksheet. You'll be taken to a completely new and fresh worksheet. We don't even have the Response Type filter anymore that we used to get rid of the preview responses. We'll need to fix that before we visualize anything.

### Applying Filters to All Worksheets

You *could*, if you wanted to, go through the same procedure of recreating the Response Type filter every time you create a new worksheet, but that's a lot of work. Luckily, Tableau gives you an easy way to apply the filter just once for all worksheets.

Go back to the worksheet where you built your "Counts by Gender" table. In the filters shelf, right click on the Response Type. From the drop-down menu, click "Apply to Worksheets" and then select "All Using this Data Source", like so:

<img src="https://tufts.box.com/shared/static/nugddkbkfatbx2dpxrse4sg0dmkpq3pk.png" alt="Applying the filter to all worksheets">

Now go to the new worksheet you just created. (It should still be called "Sheet 3" since we haven't renamed it yet.) If the Response Type filter hasn't appeared already, it will appear the moment you start using data from our data set. You can test this by dragging any one of your dimensions or measures to the Rows or Columns shelf, and watch for the filter to appear. Try it, and when you've seen it work, remove the pill you added by right-clicking on it and selecting "Remove".

 In the future, when you create new worksheets, this filter will continue to appear as soon as you start using data from this data set. 

## Charting the Responses by Gender

### Measures and Dimensions, Revisited

In the previous lessons, we placed **dimensions** (categorical variables) in the Rows and Columns shelves to create the columns and rows of a table, and we put **measures** (numerical variables) in the body of our table using the Text card in the Marks shelf.

Actually, the Row and Columns shelves are a little more versitile than that. They can handle measures too, but Tableau won't interpret measures as the columns or rows of table. Instead, it's going to assume you're trying to create a visualization. 

Try it now. Drag the "SurveyData(Count)" Measure to the Rows shelf. Tableau will notice that it's a measure and immediately try to create a bar chart:

<img src="https://tufts.box.com/shared/static/kpkgjemyzv30hejqf7xn7zgcobdmyr1c.png" alt="Bar chart with a single column">

### Using Dimensions to Create Multiple Bars in the Bar Chart

At the moment, it's a bar chart with a single column. If we want to split up our data into multiple bars, we're going to need a dimension. Let's use our Gender dimension for this purpose. Drag Gender from the  Data pane to the Columns shelf. Your chart should now look like this:

<img src="https://tufts.box.com/shared/static/j2fi6spcqid6aue8uxorso0w9dtnfpqt.png">

Congratulations, you just made your first visualization! Let's rename the worksheet now so we don't forget what we created in it. Right-click on the worksheet tab, which should currently same "Sheet 3", and rename it something like "Gender Counts Bar Chart". 

> **Changing the Visualization Type (Optional Exercise)** 
>
> There are often several different ways of visualizing the same data, and this is a good opportunity to try some of them out. 
> 1. Duplicate the worksheet so you don't lose the current version of the visualization you just made.
> 2. In the duplicated worksheet, take note of the panel of visualization options on the right side of your worksheet. Note that hovering your mouse over each visualization style gives the name of the visualization at the bottom of the panel.
> 3. We suggest you try the treemaps, the packed bubbles, and the pie chart. When you click on any of these, take note of how Tableau rearranges the Gender and CNT(SurveyData) pills to conform to the new chart type. Feel free to move things around and see what happens; you can always use Ctrl+Z to undo anything, or go back and re-duplicate the Gender Counts Bar Chart to start this exercise over. 

## Charting the Number of Respondents for Each Rating

Now let's try something slightly more tricky and see if we can create a similar-looking chart for the number of survey respondents corresponding to each overall rating that was given for Jumbo's Croissants. 

First, let's create a new worksheet and call it "Counts by Rating Chart". In the new worksheet, we'll do the same thing as before and drag the "SurveyData(Count)" measure to the Rows shelf. As before, you'll see a single-column bar chart appear, in addition to the filter that you applied earlier to all worksheets using this data set. 

Now remember that for our previous chart (the one of counts by gender), we got multiple bars by adding Gender to the Columns shelf. Let's try that again, but this time using the Overall Rating measure instead of the Gender dimension. The result should look like this:

<img src="https://tufts.box.com/shared/static/cdoz5bnahga4w3pvn0wjl0j3jicf139r.png" alt="Counts by rating chart, incorrectly done with two measure fields">

What happened? Because Overall Rating is a measure, Tableau doesn't see it as a viable way to group your data. Since you can't make a bar chart with two measures, it made a guess that you wanted to make a scatter plot instead; note that the chart type in the "Show Me" menu at the right now has the scatter plot selected. (You might also curious as to why it chose to sum up your overall rating measure when you added it to the columns shelf. See the following note for an explanation.)

> **Note:** Tableau won't *always* look for a summary statistic (e.g. a sum or average) for measures placed in the Columns or Rows shelves or on one of the Marks cards. That depends on the type visualization or chart you're trying to make, as well as on what data you've already added to these shelves. In this example, the overall rating in the Columns shelf had to work with the CNT(SurveyData) you'd already placed in the Rows shelf. Since a count is just a single data point, it had to find a way to reduce the overall rating to a single data point as well, because in a scatter plot the number of data points has to be equal along each axis.

This would all be very helpful if we were trying to create a scatter plot. We wanted to created a bar chart, however, with different columns for each rating from 1 to 10 and the number of people that responded who gave each rating. But clearly we can't use the overall rating the same way that we used gender because it's not a dimension. How do we solve this?

### Converting a Measure into a Dimension

There are a number of ways to solve our problem, but the easiest is most likely just to convert the Overall Rating into a dimension. This won't always work, but it works here because the Overall Rating only takes a finite set of values that would work well when treated as categorical. 

There's just one catch; if you did the optional lesson where you calculated average ratings by gender, then you've already used the Overall Rating field as a measure, so we don't want to mess up our previous work by converting it now. And even if you didn't do the optional lesson, it's possible you might still want to do calculations on the Overall Rating. To avoid any problems, let's create a copy of the overall rating field and make the copy into a dimension, without altering the original. 

To make a copy of a field, right the field and select "Duplicate" from the drop-down menu, like this: 

<img src="https://tufts.box.com/shared/static/x1aksh520fwf7858rsi4ciczi90ewtgf.png" alt="Duplicating the Overall Rating measure">

A new measure will appear called "Overall Rating (Copy)". Let's rename that. Right-click on the copy and select "Rename", and then call it "Overall Rating Dimension". 

Now, right-click on "Overall Rating Dimension" and select "Convert to Dimension from the drop-down menu. 

<img src="https://tufts.box.com/shared/static/033c1915b247owcp5p1u5m5qdtwt62tm.png" alt="Select Convert to Dimension from the drop-down menu">

You will now see that "Overall Rating Dimension" has been moved to the dimensions grouping within your data pane (above the horizontal line).

Now let's replace the measure version of Overall Rating in your Columns shelf with the dimension one. Right-click on SUM(Overall Rating) in Columns shelf and select "Remove", and then drag the overall rating dimension that you just created to the Columns shelf.

You should now see the chart we were looking for: 

<img src="https://tufts.box.com/shared/static/8h5n3i6zfbbbijuzwy8vru3dlx2o0el8.png" alt="Bar chart of counts by rating after Overall Rating was converted to a dimension">

> **Tip:** This worked because our Overall Ratings field worked equally well as a measure or as a dimension. A more general solution, which would work for measures that are continuous (i.e. contain a range of possible values that don't convert well into categories) would be to use a histogram instead of a bar chart. Histograms group similar measure values into "bins" which can function similar to a dimension. For more information on histograms and bin variables, see the official Tableau documentation. 


### Cleaning Up Your Bar Chart

Your bar chart is starting to look good! Let's tweak a few things.

#### Retitling your Chart

By default, the title you see at the top of your chart is whatever you named your worksheet at the tab on the bottom of your screen. Let's change it to something else. 

Right-click on the title and select "Edit Title" from the pop-up menu. 

<img src="https://tufts.box.com/shared/static/ij8it2bx6ypebai5wl4caiif7zutsd58.png" alt="Select edit title from the pop-up menu">

A window will appear with options to change the text. You will see that there's already some text written there: \<Sheet Name\>

The angle brackets (i.e. the less than and greater than symbols) are a special symbol in Tableau's markdown language to indicate that this isn't normal text. In this case, they're being used to tell Tableau to retrieve the name of the worksheet. If you want to change the text for the title, you need to delete everything here, *including* the angle brackets, and replace it with the text you want. Let's change it to "Number of Responses by Rating"

<img src="https://tufts.box.com/shared/static/7c6xj8oe7ltvgvn2rdisczjqjxboph6x.png" alt="Title text window">

Feel free to change the font or other stylistic properties of the text if you like (just make sure to highlight the text you want to change first). When you're ready, click "OK". 

#### Remove Horizontal Axis Label

Tableau has labeled the horizontal axis for us using the name of the field we used to create it—in this case, the Overall Rating Dimension field. Since it's fairly obvious that these are ratings, we can delete this and get rid of visual clutter. Right click on the words "Overall Rating Dimension" and select "Hide Field Labels for Columns" to delete it.

<img src="https://tufts.box.com/shared/static/zqv5ycu59rmy528wim3ex8uws4rgbg8c.png" alt="Remove unnessessary field heading">


#### Edit Vertical Axis Label

We can also edit the vertical axis label to look a bit more professional. 

Right-click the vertical axis where it says "Count of SurveyData", and select "Edit Axis" to get the Edit Axis pop-up window. Under Axis Titles, there's a field called "Title" that contains the text used to label the axis. Where it says "Count of SurveyData", change it to something like "Number of Responses". 

## The Completed Chart

If you've followed along until this point, you should now have something that looks like this:

<img src="https://tufts.box.com/shared/static/zok2jsv9swljavgf9ky27kc7d3kea1rk.png" alt="Completed Counts by Rating Chart">

## Next Steps

You may only want to use Tableau to create charts and tables that can be used in presentations, in which case you've already covered the main topics you'll need to know. 

Tableau, however, has another feature that allows you to bring the contents of multiple worksheets into a single, attractive, and interactive layout known as a **dashboard.** Dashboards are great when you want to present your data on a website and allow your end-user to explore your data with additional menus, mouseover text, and other interactive elements. 


## Additional Resources

- [Tableau Desktop and Web Authoring Help (Official Tableau documentation))](https://help.tableau.com/current/pro/desktop/en-us/default.htm)