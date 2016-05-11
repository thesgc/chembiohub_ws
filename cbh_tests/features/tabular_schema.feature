Feature: A tabular data schema exists on the project API with various fields to help with rendering the handsontable

    
    Scenario: Tabular Data Schema included in fields all exist in the schema part
        Given I create a project as before
        Then I can list the projects on the system
        Then the tabular schema included in fields all exist in the schema


    Scenario: Tabular Data Schema included in fields all have a knownby value in the schema part
        Given I create a project as before
        Then I can list the projects on the system
        Then all the tabular fields have a knownBy value



    Scenario: Tabular Data Schema included in fields have a project specific schema with a string ID inside it if they are custom fields
        Given I create a project as before
        Then I can list the projects on the system
        Then all the tabular Data Schema included in fields have a project specific schema with a string ID inside it if they are custom fields