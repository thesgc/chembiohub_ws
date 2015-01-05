

Chembl integration -

We have got the database working locally and can save things to it

Points of interest -

- To be more productiove, we need to make the database migration system of django work and to clear out a lot of the custom code that is now effectively not needed.
- This can come late and will involve working together with EBI to upgrade the code to Django 1.7

- Salts can have there own records and parents in Chembl, these are just not shown to the outside world in the standard chembl export
- I now understand how I will generate random string ids

Unit test

- I have demonstrated a unit test to save data to my own batch table which includes custom fields using HSTORE columns
- This will be integrated today and tomorrow so that by Thursday, at least the single molecule upload works, ideally the batch upload too

Permissions and groups

- We have not nailed down the exact mechanism for permissions and groups and we aim to make contact with more IT people this week to get group data from departments
- Karl will be providing group data from Chemistry as a CSV export with appropriate fields

Weird salts, spiro etc.
