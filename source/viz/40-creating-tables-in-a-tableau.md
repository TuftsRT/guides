# Creating a Table in a Tableau Worksheet

Tableau is a data visualization tool, but it can be good for creating summary tables as well. This is often a great first step for understanding your data, and you may even choose to include some summary tables in your final dashboards.

## Objective

- Generate a table to summarize your data in a Tableau worksheet

## Getting Started with Worksheets

If you are following along with the data used in this tutorial, you should have already loaded the survey data into Tableau. Now it's time to create some summary tables for our data, but first, we need to navigate to a new worksheet.

In Tableau, **worksheets** are places to create tables or visualizations. Each table or plot will require its own separate worksheet. Later, we will show you how to combine the plots and tables from your worksheet into a dashboard.

Creating a worksheet is easy. In fact, Tableau has already created the first one for you. You will see it at the bottom of your screen:

<img src="https://tufts.box.com/shared/static/g4bfgtzsm8yznyy29l7o090xvxvj3r7p.png" alt="Go to Sheet 1">

Tableau has helpfully highlighted the worksheet in orange and added a "Go to Worksheet" prompt to guide you. Click on "Sheet 1" to go to the worksheet. Later, when we need to create more worksheets, we will use this button:

<img src="https://tufts.box.com/shared/static/fiz7k0izjgomih92uodqzjrs2ss32ovn.png" alt = "Create new worksheet button">

## Creating Your First Table

Let's create a table! The following tutorial introduces you to each of the panels and tools that you will commonly use in a Tableau worksheet. Along the way, we will show you how to use each tool to create a table that shows summary statistics of our data by gender.

### The Worksheet Interface

When you open a new worksheet, your screen should look like this:

<img src=https://tufts.box.com/shared/static/sjzkv84ck2vnih4afikt1eeithg634vp.png alt = "New worksheet interface">

### The Data Pane

The first thing you may notice is that Tableau has listed all the fields (columns) from your data source in the **data pane**, which is located in the **side bar** on the left-hand side of your screen. You might have also noticed that there are a few new fields as well that Tableau has created to help you out. We'll cover these shortly.

Note that the fields listed in this pane are divided into two groups separated by a horizontal line:

- **Dimensions**, above the line with blue icons, correspond to categorical data that are useful for defining the rows and columns of your table.
- **Measures**, below the line with green icons, correspond to numerical data that can be summed up or mathematically manipulated to create summary stats in your tables.

### The Columns and Rows Shelves, Part I

Tables are organized in columns and rows, which can be controlled by the "Columns" and "Rows" panels near the center-top of your screen. Actually, Tableau uses slightly different terminology here. Instead of calling these "panels", Tableau calls them "shelves". There are other shelves in your worksheet too, such as the "Filters" and "Marks" shelves, which we'll get to shortly.

For now, let's try using our data to define the rows of a table. In your Data panel, locate the "Gender" field. Click and drag it to the "Row" shelf. You should now see something that looks like this:

<img src=https://tufts.box.com/shared/static/gt8dm7nbotipvahryotz3vekfm2dmllj.png alt="Table rows defined by gender">

"Gender" now appears as a colored, pill-shaped element in the Rows shelf. The shape gives these elements their name: **pills**. Pills are created when you move fields and calculations into shelves, and they can easily be dragged around between shelves to change the structure of your tables and visualizations.

### The View Panel

When you placed the Gender pill in your row shelf, you probably noticed that your table appeared in the **view**. This is the part of your worksheet interface where you get to see your tables and visualizations come to life.

You'll notice, however, that it's not much of a table yet. We have row names, but we haven't told Tableau what to put in the table itself. Tableau has filled in the body of the table with "Abc" as a placeholder while waiting for you to tell it what data to use here.

### The *SurveyData (Count)* Field

Let's add some data to the table. To start with, it would be useful to know how many observations there were for each gender.

Tableau has actually created a special field just for this purpose. Take a look at the data pane and locate the field called \*SurveyData (Count)". This field will tell us the number of observations in our data, or in any breakdown of our data defined by our table.

But where should we put it so that it goes into our table? This is where the Marks shelf comes in.

### The Marks Shelf

The "Marks" shelf is a powerful tool for controlling different attributes of your tables and visualizations in Tableau. Take a minute to locate it to the right of your data pane, just below the shelves marked "Pages" and "Filters". Notice that it contains a number of cards that can be used to adjust various aspects of your table.

#### The Text Card

At the moment, we're most interested in the "Text" card, which determines the content of our table. From your Data pane, click and drag the *SurveyData (Count)* field to the Text card. Once you've dropped it on the Text card, it will reappear as a pill next to a little text icon in your Marks shelf, as shown in the image:

<img src=https://tufts.box.com/shared/static/f5ff00tjdu2erbzy5hx65al99zo0yfv4.png alt="Observation counts added to the table" >

Now take a look at your table in the view panel. Did you notice that your table now has values in it? We now know how many respondents of each gender were in our survey data.

#### The Color and Size Cards (Optional)

If you like, you can also try playing around with the "Color" and "Size" cards. Try clicking on each and playing around with the options that pop up to change the color and size of the text in your table. If you don't like your changes, you can always undo them by pressing Ctrl+Z on a Windows PC or Command+Z on a Mac.

You can also allow attributes like color to be determined by your data. Try clicking and dragging the Gender field onto the color card to see your data color-coded by gender. It is important to consider the careful use of color in your visualizations.

#### The Tooltip Card (Optional)

Tableau is meant to produce interactive tables and visualizations. One way to do this is by having a message pop up when you hover your mouse over your table (or visualization) in the view panel.

By default, Tableau constructs this for you. If you hover your mouse over the numerical values in your table, you will already see the pop-up message appear with some information about your data. If you want to play around with the appearance or content of this pop-up message, click on the "Tooltip" card. For more information on how to modify Tooltip messages, consult the official Tableau documentation.

### The Column and Row Shelves, Part II

So far we've split up our data by gender. But why stop there? Tableau allows you to split up your data by multiple dimensions simultaneously.

#### Tables with Both Rows and Columns

Let's try adding our "Response Type" field as a dimension in our table. Actually, there are several different ways to do this. Let's try dragging it to our columns shelf first, like so:

<img src=https://tufts.box.com/shared/static/sjqtr2fdofri91jjo4i7kng2ii2bvn03.png alt="A table with both columns and rows">

As you can see, a column has been created for every value of the "Response Type" field.

#### Tables with Nested Rows

The Rows and Columns shelves can each support more than one pill. When you place multiple pills on either the Rows or Columns shelves, Tableau will *nest* them, one within the other, according to the order that you put them in.

For example, let's try dragging our Response Type pill down from the Columns shelf to the row shelf. If you put it *after* Gender, you will get a result that looks like this:

<img src=https://tufts.box.com/shared/static/4hv4hxdrias8fsn5b78nvkm4fpdfd68l.png alt="Response Type Nested in Gender">

Whereas if you place the Response Type field *before* Gender, you will see something like this:

<img src=https://tufts.box.com/shared/static/umtzgxvrtcz9rso0su42tnynzlo082cq.png alt="Gender Nested in Response Type">

In our example, since we have only two dimensions, assigning Gender to rows and Response Type to columns is probably the cleanest way to organize our table. However, you will often want to add three or more dimensions to a table, so knowing how to nest data effectively can come in handy.

### The Filter Shelf

In the above section, we showed how we could get multidimensional tables by adding fields to the Columns and Rows shelves. But do we really need to see the Response Type? The data from the Survey Preview responses aren't very useful, since these were just employees testing the survey. It might be more useful just to get rid of them entirely.

This is where the Filters shelf comes in handy. Take a second to locate it just above the Marks shelf. The Filters shelf helps us decide which rows of our data are to be included in our table and which are to be omitted.

Let's remove Response Type from our table and use it to create a filter. Right-click the Response Type field in your Columns or Rows shelf and select "Remove" to go back to the single-dimensional table made up only of gender.

Now, from your Data pane, click and drag "Response Type" to the "Filters" pane. You will now see the following menu:

<img src=https://tufts.box.com/shared/static/1jqsbq834c1hh8m7bw7o30d8s9dpcymt.png alt="Filter Menu">

These are the values of the Response Type field. We want to keep only the "QR Code" respondents, since these are the customers who saw the QR Code posted at the bakery and decided to take the survey. Make sure "QR Code" is checked and "Survey Preview" is not checked. Then press "OK".

If you look at your data table now, you will notice that the numbers in the table have gone down, reflecting the fact that you have omitted some of the rows of the original data set.

### Aliases

Chances are, you've probably guessed that "F" stands for "Female", "M" stands for "Male", and "NB" stands for "Non-Binary", but it would look a lot more professional if these words were spelled out. Aliases are a convenient way to give better display names to values in your data without having to change the dataset itself.

Let's give "F" the alias "Female". In your table, right click on "F" to reveal a drop-down menu and select "Edit Alias..."

<img src=https://tufts.box.com/shared/static/yb48ojkgbddfoepadq4ksfnsq3qtkkog.png alt="Right Click to Select Edit Alias">

A prompt will appear in which you can type in the alias. Type "Female" and press ok. Do these steps for "Male" and "Non-Binary" respondents as well, and your table should now look like this:

<img src=https://tufts.box.com/shared/static/zhqae0mmsdf79qznwixice49os0pimfp.png alt="Gender Table with Aliases">

## Finishing Up

### Rename the Worksheet

It's good practice to name each worksheet so you can keep track of its contents. To do this, go to the worksheet tab at the bottom left-had corner of your screen. Right-click the worksheet tab and select "Rename". Rename your worksheet "Counts by Gender" or whatever name feels best to you.

### Optional: Remove Row Heading

Good visualizations often benefit from visual simplicity, and sometimes it can help create a cleaner look if you remove unnecessary information. In our table, for example, it's pretty obvious that the left-hand column contains gender information, so having a column header might just be unnecessary clutter. If you like, you can hide it by right-clicking on the column heading, right where it says "Gender", like so:

<img src="https://tufts.box.com/shared/static/zu7hiezssc6x6pxew6xyh3vdxphmmij1.png" alt="Test">

You will now see a more stream-lined version of your table with the column heading removed.

## Next Steps

This lesson showed you how to create a simple summary table in Tableau. If you want to create more complex tables in Tableau, with multiple measures or more complex calculations, we highly encourage you to work through the next lesson on [creating more complex tables](<>). This lesson is optional, however, and if you prefer, you can skip straight ahead to the following lesson on [creating your first data visualizations.](<>)

## Additional Resources

- [Tableau Desktop and Web Authoring Help (Official Tableau documentation))](https://help.tableau.com/current/pro/desktop/en-us/default.htm)
