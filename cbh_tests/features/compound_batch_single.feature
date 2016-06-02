Feature: A user can save a compound with and without a structure and can clone a compound by removing the IDs and re submitting

    Scenario: A compound batch is created from a drawing
        Given I start the qcluster
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
        When I stop the qcluster
        Then I see the right qcluster output


    Scenario: A compound batch with structure has a UOx ID  in the uuid field
        Given I start the qcluster
        Given I have loaded the fixtures for project types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        When I list compound batches in the system
        then I see no compound batches
        Given A compound batch is created from a drawing as before
        When I list compound batches in the system
        Then the created compound batch has a uox id in the uuid field
        Then the created compound batch has a multiple_batch_id
        When I stop the qcluster
        Then I see the right qcluster output

    Scenario: An inventory record is created
        Given I start the qcluster
        Given I create a project as before
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add the blinded batch id to my compound POST data as EMPTY_ID
        and I add the project primary key to the compound data
        When I submit the compound by POST request
        Then a compound batch is created
        and the blinded batch ID is generated
        When I list compound batches in the system
        Then the created compound batch has a uox id in the uuid field
        Then the created compound batch has a multiple_batch_id
        When I stop the qcluster
        Then I see the right qcluster output


    Scenario: I can request a single compound batch via the get detail api
        Given I start the qcluster
        Given I create a compound batch from a drawing as before
        When I request the compound batch with ID 1 from the get_detail api
        Then the batch response is OK
        When I stop the qcluster
        Then I see the right qcluster output









