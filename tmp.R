
mock_abun <- seq(0.1, 0.3, by = 0.05 ) * 100;

mockabun <- c( mock_abun, mock_abun, 0 ) / 2;

species <- c( "Enterobacteriaceae", "Lachnospiraceae", "Bacteroidaceae", "Lactobacillaceae", "Clostridiaceae", 
              "Ruminococcaceae", "Prevotellaceae", "Erysipelotrichaceae", "Streptococcaceae", "Enterococcaceae", "Other" );

speabun <- tibble( id = character(), genus = character(), abundance = double() );

for( i in LETTERS[1:10] ){
  speabun <- bind_rows( speabun, 
                        tibble( id = i, genus = species, abundance = sample( mockabun, 11, replace = FALSE ) ) );
}

write_tsv( speabun, path =  "data/talk09/mock_species_abundance.txt" );

