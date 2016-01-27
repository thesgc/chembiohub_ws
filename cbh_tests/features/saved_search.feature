Feature: A user can save a search and owns it, organisations can also own a saved search

    Scenario: A saved search is created
        Given I have loaded the fixtures for project types and data types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding the saved search project type and the static custom field config schema for saved searches
        When I POST a project to cbh_projects
        Then the project is created
        Given A URL to redirect the user to and a GET request URL and an alias and description for my saved search
        and I add the blinded batch id as EMPTY_STRING
        and I add the project key
        When I send the search by POST request
        Then The saved search response is created

  #  Scenario: A saved search can be assigned to an organisation by any organisation member
  #      Given A search has been saved as before
   #     When A form is used on the front end to build a request with a organisation to own the saved search
   #     and I send the saved search permission by post request
   #     Then the response is updated

   # Scenario: A person from an organisation can list saved searches for that organisation
     #   Given an organisation owns a saved search as set up before
    #    and a separate organisation member logs in
      #  When an organisation member lists saved searches
      #  Then they see the existing saved search


    #Scenario: The owner of an organisation can remove a saved search from an organisation



   # Scenario: A member of an organisation cannot remove a saved search from an organisation
