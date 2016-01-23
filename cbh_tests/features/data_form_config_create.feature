Feature: A user can create a data form config for assayreg

Scenario: Add custom field configs
    Given I have loaded the fixtures for project types and data types
    Given testuser has the cbh_core_model.add_project permission
    When I log in testuser
    Then I can list the data types on the system and there are 4
    Given I POST 4 custom field configs
    When I list the custom field configs on the system there are 5
    Given I create a data form config using the 4 custom field configs
    When I post the template data form config
    Then The data form config is created