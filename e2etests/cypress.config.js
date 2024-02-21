const { defineConfig } = require("cypress");
const cucumber = require('cypress-cucumber-preprocessor').default;
module.exports = defineConfig({
  'cypress-cucumber-preprocessor': {

    nonGlobalStepDefinitions: true,

    step_definitions: './cypress/e2e/**/*.feature',

  },
  e2e: {

    baseUrl: 'https://www.stg.phase-mapper.matr.io/',
    username: '',
    password: '',
    setupNodeEvents(on, config) {
      on('file:preprocessor', cucumber());
      
    },
    testIsolation: true,
    supportFile: 'cypress/support/commands.js',
    specPattern: 'cypress/e2e/**/*.feature',
    experimentalWebKitSupport: true,
  },
});
