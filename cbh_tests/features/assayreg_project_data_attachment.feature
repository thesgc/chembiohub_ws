Feature: Files can be attached to the project data component allowing a flexible file field

    Scenario: A flowfile object is created from the front end
        Given I upload a file to django flowjs
        When I GET the file object via the URI from the (made up) session ID
        Then The file has been created and the file details can be retrieved from the server


    Scenario: An attachment object is saved by an owner of the project
        Given I create a datapoint classification as before
        and A flowfile object is created from the front end
        When I list the nested datapoint classifications in the project
        and I combine the flowfile id with the parent data point classification
        and POST the attachment object to the base attachment API
        Then the attachment object has been created
        When I make a request to the raw image url
        Then the raw url delivers a response in png format


    Scenario: The attachment object is saved by an editor of the project
        Given I create a datapoint classification as before
        and A flowfile object is created from the front end
        and I remove all of the testusers permissions
        and I make testuser an editor of the first project in the list
        When I list the nested datapoint classifications in the project
        When I combine the flowfile id with the parent data point classification
        and POST the attachment object to the base attachment API
        Then the attachment object has been created
        When I make a request to the raw image url
        Then the raw url delivers a response in png format


    Scenario: The attachment object cannot be saved by a viewer of the project
        Given I create a datapoint classification as before
        and A flowfile object is created from the front end
        and I remove all of the testusers permissions
        and I make testuser a viewer of the first project in the list
        When I list the nested datapoint classifications in the project
        When I combine the flowfile id with the parent data point classification
        and POST the attachment object to the base attachment API
        Then the attachment object is not created as unauthorized


    Scenario: The attachment object can be viewed by a viewer of the project
        Given An attachment object is saved by an owner of the project
        and I remove all of the testusers permissions
        and I make testuser a viewer of the first project in the list
        When I make a request to the raw image url
        Then the raw url delivers a response in png format


    Scenario: The attachment object cannot be viewed by someone without permissions
        Given An attachment object is saved by an owner of the project
        and I remove all of the testusers permissions
        When I make a request to the raw image url
        Then the request to the raw url is rejected as unauthorized 

    Scenario: The attachment object cannot be viewed if the user has logged out
        Given An attachment object is saved by an owner of the project
        and I log out
        When I make a request to the raw image url
        Then the request to the raw url is rejected as unauthorized 