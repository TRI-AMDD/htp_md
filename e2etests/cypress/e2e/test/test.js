import {Given, When, Then} from "cypress-cucumber-preprocessor/steps";

Given("precondition",()=>{
    cy.visit("https://www.google.com")
})

When("action",()=>{
    cy.title().should('eq','Google')
})

Then("testable outcome",()=>{
    cy.title().should('eq','Google')
})