Feature: A user can save a search and owns it, organisations can also own a saved search

    Scenario: A saved search is created
        Given A URL to redirect the user to and a GET request URL
        When I send the search by POST request
        Then The saved search response is created
        and the creator of the saved search is set as its owner

    Scenario: A saved search can be assigned to an organisation by any organisation member
        Given A search has been saved as before
        When A form is used on the front end to build a request with a organisation to own the saved search
        and I send the saved search permission by post request
        Then the response is updated

    Scenario: A person from an organisation can list saved searches for that organisation
        Given an organisation owns a saved search as set up before
        and a separate organisation member logs in
        When an organisation member lists saved searches
        Then they see the existing saved search


    Scenario: The owner of an organisation can remove a saved search from an organisation



    Scenario: A member of an organisation cannot remove a saved search from an organisation
