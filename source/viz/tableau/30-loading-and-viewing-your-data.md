# Loading and Viewing Your Data

## Objectives

- Load and view data from a Microsoft Excel spreadsheet into Tableau
- Determine and set data types (numerical or text) for imported data

## Link to Tutorial Data

For this tutorial, we will be using some sample data built to resemble the output of a hypothetical customer research survey. Suppose that Jumbo, excited by his experience at Tufts's campus in Talloires, France, has decided to start a baking business selling croissants. Eager to get feedback on his baking, he has set up an online survey for participants to rate their experience on a 1-to-10 scale and provide some basic demographic information.

If you wish to follow along with the tutorial using the data that Jumbo collected from his customers, click on the following link or paste it in your browser:

https://tufts.box.com/v/survey-data-tableau

Download and save the file to your preferred location on your computer.

## Launching Tableau and Loading Your Data

When you launch Tableau Desktop, you’ll see the **Welcome Screen**, which includes options to:

- Connect to a data source.
- Explore tutorials and sample workbooks.
- Open recent workbooks.

Here you will see options to load many different data types. If you are following along with this tutorial, choose "Microsoft Excel" and load it from wherever you saved your dataset.

<img src="https://tufts.box.com/shared/static/cztzr3k03vnb5c3zp4c0k10nw00bbdr4.png" alt="Tableau Welcome Screen" width=60%>

You'll note that Tableau supports many types of data, including spatial (GIS) files, common statistical file formats from Stata, SPSS, and SAS, and even connecting directly to servers and databases. In this walkthrough, we will use one of the most common data sources: Microsoft Excel files.

## The Data Source Screen

Once you have loaded your Microsoft Excel file, you will be directed to the Data Source screen, which should look like this:

<img src="https://tufts.box.com/shared/static/jdpdmeijrxvvlmne72podd7pevf1p008.png" alt="Data Source Screen" width=100%>

This screen was designed to allow users to upload one or more sources of data and manage the relationships between them. For this tutorial, we are working with only a single data source. Advanced users may wish to consult the official Tableau documentation for more information on managing multiple datasets.

On the left-hand side, note the "Connections" panel that lists the data you have loaded into Tableau. You should see your Excel file listed here as it is in the example.

### Connections

Below the connections panel, you will see a panel labeled "Sheets". This contains a list of all the Excel sheets Tableau identified in the file you uploaded. We need to tell Tableau which sheet we're working with. Drag and drop the Excel sheet called "SurveyData" to where it says "Drag tables here."

Tableau will take a moment to pull the information from this Excel sheet. Once it has finished, your Data Source screen should look like this:

<img src="https://tufts.box.com/shared/static/07okvpnlxflorif6y19n9cnql76pyr6m.pnghttps://tufts.box.com/shared/static/07okvpnlxflorif6y19n9cnql76pyr6m.png" alt="Data Source Screen with Worksheet Loaded" width=100%>

A new panel has appeared that displays the content of the Excel sheet. On the left side of this panel, you will see a list of "fields", or columns in your dataset. On the right, you will see a preview of the dataset itself.

Take a moment to look at the data. In the "Response Type" column, you may notice that most of the data corresponds to respondents who scanned a QR code to take the survey. The first two entries, however, are "Survey Preview" responses, which were most likely "fake" responses generated when the survey creator was testing out the survey. We will need to remember these responses and filter them out later.

### Data Types

#### Text versus Numeric

When Tableau imports your data, it tries to determine whether each column should be treated as text, numbers, or something else. The data type chosen by Tableau can be seen as a little blue icon in the "Type" column of the "Fields" table, and can also be seen just above the column names in the data preview. Text columns, also known as "strings", are denoted by "Abc", whereas numerical data are denoted a pound sign, #.

#### Why data types matter

Depending on the data type designation, Tableau will treat your data differently. For example, it may automatically try to calculate sums of numerical data when included in certain tables or visualizations. Sometimes, however, you may disagree with the column type chosen by Tableau. In our case, for example, it may make more sense to treat the Response ID as text rather than a numerical value, given that you would never need to sum up the ID numbers or do any kind of mathematical operation on an ID.

#### Changing data types

To change the data type for the ResponseID variable in our data, simply click on the data type icon as shown in the image below. This will reveal a drop-down list showing the different data types that Tableau supports. Most of these types will look familiar to you, such as numbers and dates. You may or may not have encountered a "Boolean" however, which is a type of data that can only take on the values of "True" or "False". This can be useful as an answer to a true-false question, or as an indicator variable to indicate whether a condition is met. (If Jumbo had asked if a respondent was a Tufts student, for example, he might use a Boolean variable and create a field called "IsStudent" that would be "True" for students and "False" for everyone else.) For more information on different data types and how to use them, consult the official documentation on the Tableau website.

In our case, we want to convert our ResponseID field to a String, which is the data type for text. Click on "String" from the dropdown list, and you will see the column's data type icon change to "Abc".

<img src=https://tufts.box.com/shared/static/04h33y6od2yunt100hetlt9kruhsu75u.png alt = "Changing a Column Data Type" width = 60% >

## Next Steps

You have finished loading your data. You only need to do this once per Tableau workbook for each dataset. Now it's time to create your first Tableau worksheet!

## Additional Resources

- [Tableau Desktop and Web Authoring Help (Official Tableau documentation))](https://help.tableau.com/current/pro/desktop/en-us/default.htm)
