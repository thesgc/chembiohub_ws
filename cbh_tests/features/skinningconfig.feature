Feature: There is a skinning config available with appropriate fields

    Scenario: App can be initialised by calls to backend
        Given I have loaded the fixtures for project types
        Given testuser
        When I log in testuser   
        Then I can see the skinning configuration
        Then all required fields are in the skinning configuration