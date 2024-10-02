library(tidyverse)

columns <- c("i", "area", "length",
             paste0("bodyAngles", 1:14), paste0("curvature", 1:13),
             "x", "y", "xHead", "yHead", "distToFoodHead",
             "distToFoodBody", "DMPevent", "noseLift", "pumpingRate", "comSpeed1",
             "comSpeed10", "bodyAxisSpeed1", "bodyAxisSpeed10", "omegaTurn",
             "autofocusing", "laserOn")

dir <- "outs"

exp_names <- list.files(dir)

dat <- tibble(exp_name = exp_names) |>
  separate_wider_regex(exp_name,
                       patterns = c(genotype = "^[zdelta]+", "_",
                                    replicate_nb = "[0-9]+", "_fullDataTable.csv")) |>
  mutate(data = map(exp_names,
                    \(exp_name){
                      read_tsv(file.path(dir, exp_name),
                               col_names = columns,
                               skip = 1)
                    })) |>
  unnest(data)



medians <- dat |>
  summarize(across(everything(), ~median(.x, na.rm = TRUE)),
            .by = c("genotype", "replicate_nb"))

medians |> select(-(3:38))

# plot everything
medians |>
  pivot_longer(-c(genotype, replicate_nb),
               names_to = "measure",
               values_to = "value") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~measure, scales = "free_y") +
  geom_point(aes(x = genotype, y = value))


# plot with time
dat |>
  filter(i < 100) |>
  ggplot() +
  theme_classic() +
  geom_line(aes(x = i, y = bodyAngles1,
                color = genotype, group = interaction(genotype, replicate_nb)),
            alpha = .5, linewidth = 1)


# look at variation
median_diffs <- dat |>
  summarize(across(everything(), ~median(abs(diff(.x)), na.rm = TRUE)),
            .by = c("genotype", "replicate_nb"))

median_diffs |>
  pivot_longer(-c(genotype, replicate_nb),
               names_to = "measure",
               values_to = "value") |>
  ggplot() +
  theme_classic() +
  facet_wrap(~measure, scales = "free_y") +
  geom_point(aes(x = genotype, y = value))

