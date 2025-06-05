# NCDSim

NCDSim är en simuleringsmodell för att göra framskrivningar av cancer och hjärt- och kärlsjukdomar i den svenska befolkningen, under olika antaganden om förekomsten av riskfaktorer i befolkningen, till exempel andelen rökare. NCDSim har utvecklats av Folkhälsomyndigheten i samarbete med Cancerfonden och Hjärt-Lungfonden.

Modellpopulationen baseras på SCB:s befolkningsregisterdata (LISA). Data om sjukdomsfall har hämtats från Cancerregistret och patientregistret hos Socialstyrelsen. Prevalenser för icke-kost-relaterade riskfaktorer baseras på skattningar från Folkhälsomyndighetens surveyundersökning Nationella folkhälsoenkäten. Kostrelaterade riskfaktorer baseras på skattningar från Livsmedelsverkets surveyundersökning Riksmaten. Relativa risker kopplade till riskfaktorerna baseras på Freisling m.fl. (2020), Garcia m.fl. (2020), Yi-Jie Wang m.fl. (2020) och Global Burden of Disease 2019 (2020).

En teknisk dokumentation av modellen finns [här](https://www.folkhalsomyndigheten.se/contentassets/4a8b81da030747a294c0491736494dac/ncdsim-simuleringsmodell-framtida-utvecklingen-icke-smittsamma-sjukdomar.pdf). Frågor om modellen kan skickas till analysenheten\@folkhalsomyndigheten.se

Modellen kan köras antingen från ett Shiny-gränssnitt eller direkt i R, se under respektive rubrik nedan. Modellen kan anpassas dels med olika argument till funktionen simulate_model(), och dels med olika parametrar som läses in från fil (se Parametrar) nedan. När modellen körs från gränssnittet kan argumenten och vissa av parametrarna justeras därifrån. En fördel med att köra modellen från gränssnittet är att man också kan visualisera resultaten där.

NCDSim utgår alltid ifrån en baslinjekörning, vilket utgör ett kalibrerat basscenario som eventuella alternativscenarier utgår ifrån. I Shiny-gränssnittet börjar man med att ladda eller simulera fram en ny baslinje genom fliken Baslinje, varefter (kompatibla) scenarier kan laddas eller simuleras från fliken Scenarier.

## Kör modellen från Shiny-gränssnittet

### Från R Studio
Öppna app.R och tryck på knappen "Run app" i menyraden ovanför skriptet. Därifrån kan gränssnittet antingen köras direkt i RStudio, eller genom en webbläsare.

### Genom R

Se till att du står i NCDSim-katalogen (ändra med setwd()) och kör därefter:

``` r
library(shiny)
runApp()
```
Gränssnittet kommer nu att öppnas i din standardwebbläsare.

## Kör modellen direkt i R (R prompt)

Modellen körs med funktionen simulate_model() som definieras i huvudskriptet ncdsim.r, som kan läsas in med source(). Funktionen returnerar simulerade utdata i form av ett data.table-objekt (en data frame med särskilda attribut, se R-paketet [data.table](https://cran.r-project.org/web/packages/data.table/index.html)). Om argumentet write_data_to_file sätts till TRUE sparas utdata även som en csv-fil i underkatalogen /Output, med namnet output_NCDSim\_timestamp.csv. Observera att utdata från baslinjen alltid också skrivs till en katalog som specificeras av argumentet baseline_parameters, och att alternativscenarier använder dessa data.

### Baslinjekörning

För att köra en baslinjesimulering, sätt argumentet is_baseline till TRUE, när simulate_model() anropas. Övriga viktiga argument är 

-   startyear: första simuleringsår (\<=2022)

-   endyear: sista simuleringsår (\<=2119)

-   baseline_parameters: sökväg till json-fil med baslinjeparametrar. Om detta argument utelämnas används standardvärdena i baseline.json i katalogen Input/defaults. 

Exempel på baslinjesimulering från 2021 till och med 2060:

``` r
source("ncdsim.R")
output <- simulate_model(startyear = 2021, 
                         endyear = 2060, 
                         baseline_parameters = 'path/to/my_parameters.json', 
                         is_baseline = TRUE)
```

### Scenariokörning

För att köra ett alternativscenario, sätt argumentet is_baseline till FALSE, när simulate_model() anropas, och ange en sökväg till en befintlig baslinjekörning med argumentet baseline_parameters. För att ändra riskfaktorsprevalenserna och hur förändringarna fasas in , använd följande argument:

-   cfact: en vektor med skalningsfaktorer för icke-kostrelaterade riskfaktorprevalenser. Utgångsvärdena är 1, och en prevalens kan t.ex. dubbleras om värdet sätts till 2.

-   cfact_food: som cfact för kostfaktorerna

-   cfact_startyear: startår för infasning av ändrade riskfaktorprevalenser

-   cfact_endyear: slutår för infasning av ändrade riskfaktorprevalenser

Här är ett exempel på ett scenario där prevalensen för obesitas sänks till hälften av utgångsvärdet och interventionen börjar 2025 och slutar år 2060.

``` r
source("ncdsim.R")

output <- simulate_model(baseyear = 2021, 
                         endyear = 2060, 
                         baseline_parameters = 'path/to/my_parameters.json', 
                         cfact = 
                          c(cfact_smoking = 1.0,
                            cfact_alcohol = 1.0,
                            cfact_inactivity = 1.0,
                            cfact_obesity = 0.5),
                         cfact_food = 
                          c(fruit  = 1.0,
                            wholegrains  = 1.0,
                            greens  = 1.0,
                            meat   = 1.0, 
                            salt = 1.0),
                         cfact_startyear = 2025,
                         cfact_endyear = 2040,
                         is_baseline = FALSE)
```


## Parametrar

Modellen läser in parametrar för relativa risker, kostnader, tillskrivningsfaktorer m.m. från en json-fil. Standardvärden för dessa parametrar är definierade i filen baseline.json i katalogen Input/defaults.


## Validering

Skriptet validation.R används för att skapa en uppsättning diagram för att validera resultaten från en simuleringskörning. Bland annat visas utdata från simuleringen jämsides med faktiska data över sjukdomsfall från Socialstyrelsen och med SCB:s demografiska prognos. Diagrammen sparas i en pdf-fil i /Output. Argumentet timestamp ska matcha ändelsen för csv-filen med utdata från simuleringen ifråga, under /Output.

``` r
source("validation.R")

validate_ncdsim(
  projectroot = "C:/path/to/ncdsim/folder",
  timestamp = "2025_02_19_20_04_3976") 
```