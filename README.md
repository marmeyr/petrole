# Projet Pétrole

Ce projet en R a été réaliser dans le cadre du cours "Économétrie des séries temporelles" de M1 IREF. L'objectif était d'effectuer des prévisions de la production annuelle de pétrole non raffiné au Canada pour les années 2022 à 2025. Pour ce faire, nous utilisons des données historiques de 1971 à 2021 pour un pays sélectionné aléatoirement parmi ceux disponibles dans un tableau. Les données sont extraites du site de l'OCDE [OCDE - Production de pétrole brut](https://data.oecd.org/energy/crude-oil-production.htm#indicator-chart).

## Approche méthodologique

Dans un premier temps, il a fallut déterminer si la série étuidée est issue d'un Processus Générateur de Données (PGD) stationnaire, d'un processus déterministe (DS) ou d'un processus stochastique (TS) en effectuant des racines unitaires. Cette étape était nécessaire pour identifier le nombre de différentiation requis afin d'obtenir une série stationnaire.

Dans un deuxième temps, nous avons procédé à la détermination du meilleur modèle à l'aide de l'analyse de l'Auto-corrélation Étendue (EACF) et des tests post-prévision.

Enfin, j'ai réalisé les prédictions et ai itéré le choix du modèle et des prédictions en incluant la valeur supplémentaire de 2022.

## Accès au code et aux explications

Le code R utilisé pour ce projet est disponible dans le répertoire du dépôt GitHub. Vous pouvez accéder au code source et l'utiliser pour reproduire les analyses.

De plus, les explications détaillées des résultats obtenus sont fournies dans un ensemble de diapositives au format PDF. Ces diapositives fournissent un aperçu clair des différentes étapes de l'analyse, des résultats obtenus et de leur interprétation.

N'hésitez pas à explorer le code et les diapositives pour mieux comprendre les détails du projet et les conclusions tirées de l'analyse.
