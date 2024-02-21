Cypress.Commands.add('navigatePhaseMapper', () => {
    cy.visit(Cypress.config().baseUrl)
    cy.log('NAVIGATING TO PHASE MAPPER WEB')
    cy.wait(10000)
    cy.deleteDownloadsFolder();
    cy.get(".login_login__Vxx-d").should('be.visible').click()
    cy.get("#signInFormUsername").type(Cypress.config().username,{force: true})
    cy.get("#signInFormPassword").type(Cypress.config().password,{force : true})
    cy.get('[name="signInSubmitButton"]').eq(1).click()
    cy.wait(3000)
})