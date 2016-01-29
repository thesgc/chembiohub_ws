Feature: A user can save a compound with and without a structure and can clone a compound by remoiving the IDs and re submitting

    Scenario: A compound batch is created from a drawing
        Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add a valid molfile to my compound data and call it sketch
        and I add the project key to the compound data
        and I set the type of the request to sketch
        and I set the state to validate
        When I submit the compound to POST validate drawn
        then the response from post validate drawn is accepted
        when I take the response from post validate drawn and post it to multi batch save
        then the response from multi batch save is created


    Scenario: A compound batch with structure has a UOx ID  in the chemblId field
        Given I have loaded the fixtures for project types and data types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        When I list compound batches in the system with get_list_elasticsearch
        then I see no compound batches
        Given A compound batch is created from a drawing as before
        When I list compound batches in the system with get_list_elasticsearch
        Then the created compound batch has a uox id in the chemblId field
        Then the created compound batch has a multipleBatchId


    Scenario: An inventory record is created
        Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add the blinded batch id to my compound POST data as EMPTY_STRING
        and I add the project key to the compound data
        When I submit the compound by POST request
        and I reindex the compound
        Then a compound batch is created
        and the blinded batch ID is generated
        When I list compound batches in the system with get_list_elasticsearch
        Then the created compound batch has a multipleBatchId



