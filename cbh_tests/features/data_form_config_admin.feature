Feature: The data form config admin creates a hierarchy of data form config objects

    Scenario: A data form config, when saved in the admin interface is turned into a hierarchy
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Given I create new custom field configs and data form configs based on the data given
        When I save the data form config via the admin ui
        Then All ancestor data form configs are present


        

        
