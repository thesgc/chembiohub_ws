Feature: CBH Compound Batch
    In order to keep track of compounds I would like to create 
    a batch of a new compound in the project I am an admin of.
    As a project stakeholder I do not want batches to be logged in 
    projects that the user does not have editor access to
#behave && cat typescript | aha > test.html

    Scenario Outline: Batches project privileges saving and validating
        Given a User
        and my user is member of a group
        and I have a valid molfile
        and a valid project exists proja
        and I automatically have editor permissions as creator
        and I remove my permissions
        and I have<editor> given <me_or_group> editor privileges for proja
        and I have<viewer> given <me_or_group> viewer privileges for proja        
        When <dologin> log in
        Then I <action> my cbh_compound_batch to proja and the <responsecode>


        Examples: Validation api
        |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
        |myself|      I         |           |           |  response code will be 202             | validate |
        |myself|       I        |           |   nt     |  response code will be 202             | validate |
        |myself|        I       |   nt     |           |  response code will be 401             | validate |
        |myself|         I      |   nt     |   nt     |  response code will be 401             | validate |
        |myself|   I do not      |           |           |  response code will be 401             | validate |
        
        Examples: Create api
             |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
         |myself|      I         |           |           |  response code will be 201            | create |
        |myself|       I        |           |   nt     |  response code will be 201             | create |
        |myself|        I       |   nt     |           |  response code will be 401             | create |
        |myself|         I      |   nt     |   nt     |  response code will be 401             | create |
        |myself|   I do not      |           |           |  response code will be 401             | create |

        Examples: Validation api group perms
        |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
        |mygroup|      I         |           |           |  response code will be 202             | validate |
        |mygroup|       I        |           |   nt     |  response code will be 202             | validate |
        |mygroup|        I       |   nt     |           |  response code will be 401             | validate |
        |mygroup|         I      |   nt     |   nt     |  response code will be 401             | validate |
        |mygroup|   I do not      |           |           |  response code will be 401             | validate |
        
        Examples: Create api group perms
             |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
         |mygroup|      I         |           |           |  response code will be 201            | create |
        |mygroup|       I        |           |   nt     |  response code will be 201             | create |
        |mygroup|        I       |   nt     |           |  response code will be 401             | create |
        |mygroup|         I      |   nt     |   nt     |  response code will be 401             | create |
        |mygroup|   I do not      |           |           |  response code will be 401             | create |



    Scenario Outline: Batches project privileges get list
        Given a User
        and my user is member of a group
        and I have a valid molfile
        and a valid project exists proja
        and I automatically have editor permissions as creator
        and a single batch exists in proja
        and I remove my permissions
        and I have<editor> given <me_or_group> editor privileges for proja
        and I have<viewer> given <me_or_group> viewer privileges for proja        
        When <dologin> log in
        Then I <action> my cbh_compound_batch to proja and the <responsecode>


        
        Examples: Get List api  - Note the viewer can list
        |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
         |myself|      I         |           |           |   response code will be 1_memberlist      | list |
        |myself|       I        |           |   nt     |   response code will be 1_memberlist      | list |
        |myself|        I       |   nt     |           |   response code will be 1_memberlist   | list |   
        |myself|         I      |   nt     |   nt     |   response code will be 0_memberlist      | list |

        
        Examples: Get api  - Note the viewer can list
        |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
         |myself|      I         |           |           |   response code will be 200     | get |
        |myself|       I        |           |   nt     |   response code will be 200      | get |
        |myself|        I       |   nt     |           |   response code will be 200   | get |   
        |myself|         I      |   nt     |   nt     |   response code will be 401      | get |



        Examples: Get List api  - Note the viewer can list group perms
        |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
         |mygroup|      I         |           |           |   response code will be 1_memberlist      | list |
        |mygroup|       I        |           |   nt     |   response code will be 1_memberlist      | list |
        |mygroup|        I       |   nt     |           |   response code will be 1_memberlist   | list |   
        |mygroup|         I      |   nt     |   nt     |   response code will be 0_memberlist      | list |

        
        Examples: Get api  - Note the viewer can list group perms
        |me_or_group|   dologin    |   editor  |   viewer  |   responsecode   | action |
         |mygroup|      I         |           |           |   response code will be 200     | get |
        |mygroup|       I        |           |   nt     |   response code will be 200      | get |
        |mygroup|        I       |   nt     |           |   response code will be 200   | get |   
        |mygroup|         I      |   nt     |   nt     |   response code will be 401      | get |








    Scenario: User submits a substance to project they have editor rights to, but the substance is already preregistered to other private and projects they do not have editor rights to. In one of the private projects the substance is marked as public.
        Given a User
        and I have a valid molfile
        and a valid project exists proja
        and a valid project exists projb
        and a valid project exists projc
        and a valid project exists projd
        and this substance is not pre-registered to proja 
        and this substance is pre-registered to projb as private
        and this substance is pre-registered to projc as private
        and this substance is pre-registered to projd as public
        and I have editor rights for proja
        and I have viewer rights for projb
        and I have no rights for projc
        and I have no rights for projd
        When I log in 
        and I validate a cbh_compound_batch to proja        
        Then the response will contain an id for a previously registered substance in projb
        and the response will contain an id for a previously registered substance in projd
        and the response will not contain an id for a previously registered substance in projc

'''
    Scenario: User submits ID of substance to make a batch, perhaps that they have received by email, but is not entitled to this substance.
        Given a User
        When I log in 
        and I have a valid molfile
        and a valid project exists proja
        and a valid project exists projb
        and this substance is not pre-registered to proja 
        and this substance is pre-registered to projb as private
        and I have editor rights for proja
        and I have no rights for projb
        When I register substance as a batch to proja with emailed ID
        Then the response code will be 400
        the response error will be Substance ID not found or you have no privileges for this ID
       

    



    Scenario: User registers a batch with a substance ID that does not exist 
        Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        and I have a substance ID falseid that does not exist
        When I create batch with falseid       
        the response error will be Substance ID not found or you have no privileges for this ID




    Scenario: User registers a batch with a substance ID that does not conform to the system pattern
        Given a User
        When I log in 
        and I have a substance ID blaah that is not valid
        When I create batch with blahh       
        Then the response will be error_invalid_ID_received

'''










    """
    Created using
    arrays = [("","do not"), ("","not"), ("","not"), ("401",)]
    rows = list(itertools.product(*arrays))
    rows = ["|\t" +  "\t\t|\t".join(row) + "\t\t|" for row in   list(itertools.product(*arrays))]
    print "\n".join(rows)
    """





    #Scenario Outline: Other project privileges, ivladifd project key

