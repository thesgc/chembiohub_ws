
Feature: Upload contains a salt
'''
    I wish to submit a compound or compounds which contain a salt and for that salt to be correctly registered.

    Scenario: User submits a substance containing a salt to proja. 
        Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid salted compound
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my cbh_compound_batch to proja
        and I create my cbh_compound_batch to proja
        and retain its salt

    Scenario: User submits a Chemdraw cdx file containing a salt to proja. 
        Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid salted compounds within a ChemDraw cdx file
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my cbh_compound_batch to proja
        and I create my cbh_compound_batch to proja
        and retain its salt

    Scenario: User submits a Chemdraw cdxml file containing a salt to proja. 
        Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid salted compounds within a ChemDraw cdxml file
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my cbh_compound_batch to proja
        and I create my cbh_compound_batch to proja
        and retain its salt

    Scenario: User submits a Chemdraw sdf file containing a salt to proja. 
        Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid salted compounds within a ChemDraw sdf file
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my cbh_compound_batch to proja
        and I create my cbh_compound_batch to proja
        and retain its salt


