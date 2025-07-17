# Creating More Complex Tables (Optional Lesson)

In the previous lesson, we created a simple table that gave response counts by gender. If you want to learn how to add additional measures in separate columns to your table, or to produce tables with more complex summary statistics such as percentages and averages, you will need to learn just a few more skills. 

This lesson is optional and may be skipped without disrupting your ability to follow along with later lesson.

## Objective

- Add calculations to our tables such as percentages and averages
- Use the "Measure Names" and "Measure Values" pills to add multiple measures to the same table

## Duplicating the Worksheet

In the last lesson, we created a simple table that counted the number of observations by gender. When we have a version of a table that we like, it can often be useful to keep a copy of that version before we try to build off it. Let's therefore create a duplicate of our table in a new worksheet and use the duplicate copy as the basis for creating our new, more complex table. This can help to quickly iterate through versions and test out ideas for visualizations and tables.

To duplicate the worksheet:

- Right-click on the worksheet created in the previous lesson ("Counts by Gender") in the Tableau workbook's tab bar.
- Select **Duplicate** from the context menu. This creates a copy of your table that we can modify without affecting the original.

<img src=https://tufts.box.com/shared/static/t6xeu2ff4r7rib6uclsocywji7sjmwjq.png alt="Duplicating the Worksheet">

You should now be looking at a duplicated worksheet with the name "Counts by Gender (2)". It would be helpful to rename this new worksheet. Right-click on the tab for the new worksheet, select "Rename" from the list of options, and type in "Gender Summary Table".  

## Adding More Measures

If you've followed along so far, you've already seen tables with multiple rows and columns. Adding a **dimension** (categorical data that can be used to break up our data into groups of interest) to the Rows or Columns shelves creates as many rows or columns as there are values for that dimension. 

But what if you wanted to add more **measures** (numerical data) to the body of our table? Suppose, for example, that we also wanted to include an average of the overall ratings supplied by our survey takers, broken down by gender?

Let's try it. Locate the Overall Rating field in your data pane and drag it over to the Text card in your Marks shelf. Your workspace should now look like this:

<img src=https://tufts.box.com/shared/static/yiaurr7vfwsbqeli4y9hc8egtobh5pxe.png alt="Gender Summary Table after dragging Overall Rating to the Text card">

When you dragged the Overall Rating field to the Text card, it created a new pill in our Marks shelf called "SUM(Overall Rating)" with a little "T" icon next to it to let us know it's associated with the Text card. Summing up customer ratings obviously isn't very useful; Tableau didn't know we wanted averages, so it made its best guess as to what aggregate statistic to use to summarize our data by gender. We'll need to fix that.

**If you look at the table, you'll notice one additional problem:** Since we now have two different data points placed on the Text card, Tableau thinks we want them both crammed into the same cells in our table. That makes the table visually cluttered and hard to read, so we will also need a way to tell Tableau to put them in separate columns with their own column headings. We'll come back to that in a moment in the section titled "moving measures to their own columns".

Let's start with the easier problem of changing the Overall Ratings sums to averages.

### Choosing the Right Aggregate Statistic for a Measure

Locate the "SUM(Overall Rating)" pill in your Marks shelf, and right-click on it. From the drop-down menu that appears, click on "Measure (Sum)" and select "Average" from the sub menu.

<img src=https://tufts.box.com/shared/static/pijojrme5rz9f7emmdb89z7lhhgguubi.png alt="Choosing the right aggregate statistic for our Overall Rating field">

Now take a look at your table. The sums should now be replaced by averages. Easy!


### Moving Measures to their Own Columns

Now for the slightly trickier problem: we want to move the respondent counts and the average ratings to their own columns with their own column headings. 

This is where "Measure Names" and "Measure Values" come in handy. You may have already noticed these in your Data pane, but if not, take a moment to find them now in the list of fields that Tableau pulled out of our data source.

The next few steps are a bit tricky to follow, so feel free to watch us do it first in the image below. The explanation and a step-by-step guide will follow immediately after.

<img src=https://tufts.box.com/shared/static/dtmey30m8hyjjve41jjdui5qlsjbw6iv.gif alt="Demonstration of using Measure Names and Measure Values">

#### Measure Names and Measure Values: Key Concepts 

First, let's explain the key concepts. Then we'll follow with a step-by-step guide.

- When you drag "Measure Values" to the Text control in your Marks shelf, it will create a new shelf, also called "Measure Values". All the data you want to put in the columns.
   - Note: once you've dragged "Measure Values" to the Text control, you can get rid of everything else that's on the text control.
- When you put Measure names in the columns shelf, it tells Tableau to create a column for every pill you placed in the Measure Values shelf.

#### Step-by-Step Guide

First set up Measure Names and Measure Values:

1. Drag Measure Names to the Columns Shelf.
2. Drag Measure Values to the Text control in the Marks shelf. 

Notice that this creates a new shelf called Measure Values just below the Marks shelf. Tableau likes to guess what you're going to do next, so it goes ahead and adds a few pills to this new shelf for you. For example, it went ahead and added response counts for you, much like the one we created in the previous lesson. It also tried to add the respondent's overall rating again, once again guessing that we wanted the sum. The average rating pill we created is still there in the Marks shelf where we put it though.  

The next few steps are about (a) ensuring that all the measures we want to keep are in the Measure Values shelf, and (b) getting rid of the extra measures that we don't need.

3. We definitely want the average rating variable we created earlier. Locate "AVG(Overall Rating)" in the Marks shelf and move it to the Measure Values shelf.
4. We don't need the "SUM(Overall Rating)" that Tableau created for us in the Measure Values shelf. Right click on that and select "Remove" to delete it.
5. Since Tableau created another count of survey respondent counts in out Measure Values shelf, we don't actually need the one we created earlier any more. Let's keep the new one in the Measure Values shelf and discard the old one in the Marks shelf. In the Marks shelf, right-click "CNT(SurveyData)" and select "Remove" to delete it.

Let's take stock of where we are now. If you did it correctly:
- In your Marks shelf, there should only be one pill, the Measure Values pill. 
- In your Measure Values shelf, there should be two pills, one for survey counts and one for average ratings.
- In your table, there should be two columns, one for each pill in your Measure Values shelf. 

Your workspace and table should now look like this:

<img src="https://tufts.box.com/shared/static/ar2gmhfhp204amfxboefoze97orteslb.png" alt="Gender Summary Table with average overall ratings">



> Tip: The Measure Values shelf will disappear if you have fewer than two pills on it, because Tableau assumes that it's no longer needed. To keep it from disappearing, drag all the pills that you *do* want to include to the Measure Values shelf *first*, and *then* remove all the pills that you don't want.

## Adding Table Calculations

Earlier, we created counts of our data by gender. But what if we also want to know what *percent* of respondents are male, female, or non-binary? 

Actually, percentages are slightly more complicated than the aggregate statistics we've seen so far, because Tableau sees this as an *additional* calculation to be done to data that's *already* been aggregated in some way. For Tableau, getting counts by gender identity came first, and then converting them into percents is an extra step for data that's already organized into a table in some way. To put it another way, a Table calculation is a way of doing additional work on data that's already in a table.

Right-click on CNT(SurveyData) where it appears in your Measure Values shelf, and then click "Quick Table Calculation" to get a list of commonly used Table Calculations. Then click "Percent of Total".

<img src=https://tufts.box.com/shared/static/5uel0q9dcn643w7dfp7mjy3bpfz9gdup.png alt="Quick table calculations menu">

You should now see that the count of survey respondents has been replaced in the table by percentages. If you look at the pill you just modified in the Measure Values shelf, you'll also notice that a little triangle icon has been added to the right side of the pill. That's how Tableau tells you a **Table Calculation** has been added. 

If we want the counts *and* the percentages, you can simply go to your data pane and drag SurveyData (Count) to your Measure Values shelf again. Now you'll have both, and your table should look like this: 

<img src=https://tufts.box.com/shared/static/a3pg61al1d2jqlz6mxg600t3lo22jn8n.png alt="Gender summary table with three columns before cleanup.">

## Final Cleanup

Let's clean up our table to make it look more professional. 

### Fixing the number of decimal places

Specifying the number of decimal places for each column is slightly more complicated than before. When you right-click on a data point inside your table, it will open the "Format Font" pane on the left-hand side of the screen. To format numbers for specific columns, you'll need to located the drop-down menu within this pane that says "Fields" and select the field you want to edit. For example, to select the "Percent of Total Count" column, follow the example in the image below:

<img src="https://tufts.box.com/shared/static/12oywowc6fmvg11ly777sjym7uu9cuox.png" alt="Selecting a column to reformat">

You'll then see a drop-down menu labeled "Numbers". Click on this, make sure that "Percentage" is selected, and then choose the number of decimal places you'd like. In the picture below, we are setting it to one, but you could also set it to zero if you want the percentages to be rounded to the nearest whole number. 

<img src="https://tufts.box.com/shared/static/vo8w4wp4vl669jukevgisef9yebs4c4s.png" alt="Formatting the percent of total count column">

The percentages for that column should now be rounded to whatever number of decimal places you chose. 

Let's go back to the "Fields" menu in the Format pane, and this time we'll select "CNT(SurveyData)". Once again, click on the drop-down menu labeled "Numbers", but this time select "Number (Custom)". This will allow you to display this column as a number instead of a percentage. From the customization menu that pops up, you can once again select the number of decimal places you prefer: 

<img src="https://tufts.box.com/shared/static/ivdayev3blnphnf5f7rc5e56bgg6rk6m.png" alt="Formatting the gender counts">

Feel free to edit the number of decimal places in the average rating column as well. 

Your table should now be looking a lot better, without excess decimal places!

### Changing the column names

Column names have "aliases" just like we saw when edited the values of the gender field to be more readable. To change a column's name, simply right-click on the column header and select "Edit Alias". A window will pop up asking you for the new column name. Let's try it with the column that's currently labeled "% of Total Count of SurveyData along Table (Down)", and rename it something more readable, like "Percent of Respondents". Feel free to try this with the other columns as well, coming up with your own names that will be more descriptive for the viewer. 

### Changing Column Widths

As a final step, let's make the columns wider as well. To change the column width, align your cursor with the edge of the column header until the cursor changes to a double-sided arrow that points both left and right. Click and drag the column to the desired width.

**Congrats!** If you have followed along up until now, your table should look nicely formatted and professional:

<img src="https://tufts.box.com/shared/static/xdmrvmbjeh3k1xjq9l0oxnpa88awu3jq.png">

There are many more ways to tweak and edit your table if you like. Feel free to try editing the color, size, and alignment of the measure values in your table by playing with the Color, Size, and Text cards in your Marks shelf, for example, or see the official Tableau documentation for further information on all the ways you can customize your tables. 

## Next Steps

Now that you've mastered the basics of tables, it's time to move on to your first visualizations. We'll be making some simple bar charts to demonstrate the basic elements of building a chart in Tableau.


## Additional Resources

- [Tableau Desktop and Web Authoring Help (Official Tableau documentation))](https://help.tableau.com/current/pro/desktop/en-us/default.htm)