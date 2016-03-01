Feature: Demonstrate how to use the mechanize browser to do useful things.
 
  Scenario: Logging in to our new Django site
 
    Given a user
    When I log in
    Then I see a json response with only the logged in users name in it


 
  Scenario: Trying to get user info without loggin in
 
    Given a user
    When I do not log in
    Then I see a 401 error

