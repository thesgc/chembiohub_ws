Feature: I can log in and the backend delivers the right data to start the app

    Scenario: User can log in 
        Given testuser
        When I log in testuser


    Scenario: App can be initialised by calls to backend
        Given testuser
        When I log in testuser   
        Then I can list the users on the system
        Then I can see the skinning configuration
