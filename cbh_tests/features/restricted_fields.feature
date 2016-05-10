Feature: Restricted Fields

        @wip
        Scenario: Data removed from output when field restricted
            Given superuser
            And testuser
            And I create a project as before
            And I add a restricted field
            And some compounds exist in the project with this restricted field
            When I give testuser read permissions on the project
            Then the output data from search does not contain the restricted field

        @wipnot
        Scenario: Two projects with same field one of them restricted to the testuser
            Given 2 projects exist 
            And 1 has a restricted field 
            And the other has an unrestricted field of the same type
            And there are compounds in both projects with values in that field
            When I give a user viewer permissions to both
            Then Users can see all data from both
            And if they search against the part restricted field
            Then they only see the project with the unrestricteed field


        @wipnot
        Scenario: Two projects with same field one of them restricted to the testuser
            Given 2 projects exist 
            And 1 has a restricted field 
            And the other has an unrestricted field of the same type
            And there are compounds in both projects with values in that field
            When I give a user viewer permissions to both
            Then Users can see all data from both
            And if they look at the pick from list from the part restricted field
            Then they only see the pick from list values from the project with the unrestricteed field


