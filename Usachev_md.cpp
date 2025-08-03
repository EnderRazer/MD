#include <fstream>
#include <iostream>
//#include <omp.h>

using namespace std;

#include <nlohmann/json.hpp>
//  Классы
#include "backups/BackupManager.h"   //Бекапы
#include "barostats/Berendsen.h"     //Баростат
#include "classes/EnsembleManager.h" //Менеджер ансамблей
#include "classes/ThreadPool.h"      //Пул тредов
#include "classes/Timer.h"           //Таймер
#include "classes/CellList.h"
#include "core/MDAlgorithms.h"       //Алгоритмы МД
#include "core/Settings.h"           //Класс настроек
#include "core/System.h"             //Класс системы
#include "generators/Generator.h" //Генератор координат и частиц (для инициализации)
#include "output/OutputManager.h"  // Вывод в файлы
#include "potentials/EAM.h"        //Потенциал погруженного атома (ППА)
#include "potentials/LJ.h"         //Потенциал Леннарда-Джонса
#include "thermostats/Berendsen.h" //Термостат Берендсена
#include "thermostats/Langevin.h"  //Термостат Ланжевена

using namespace std;

// Фабрика потенциала
std::unique_ptr<Potential> createPotential(const json &config) {
  std::string type = config.value("type", "LJ");
  if (type == "LJ") {
    json pot_params = config["params"];
    return std::make_unique<LJ>(pot_params);
  } else if (type == "EAM") {
    json pot_params = config["params"];
    return std::make_unique<EAM>(pot_params);
  } else {
    throw std::runtime_error("Неизвестный тип потенциала: " + type);
  }
}

// Фабрика термостата
std::unique_ptr<Thermostat> createThermostat(const json &config,
                                             Settings &settings, int pn) {
  std::string type = config.value("type", "BRNDSN");
  if (type == "BRNDSN") {
    json thermostat_cfg = config;
    return std::make_unique<ThermostatBerendsen>(thermostat_cfg, settings);
  } else if (type == "LNGVN") {
    json thermostat_cfg = config;
    return std::make_unique<ThermostatLangevine>(thermostat_cfg, settings, pn);
  } else {
    throw std::runtime_error("Неизвестный тип термостата: " + type);
  }
}

// Фабрика баростата
std::unique_ptr<Barostat> createBarostat(const json &config,
                                         Settings &settings) {
  std::string type = config.value("type", "BRNDSN");
  if (type == "BRNDSN") {
    json barostat_cfg = config;
    return std::make_unique<BarostatBerendsen>(barostat_cfg, settings);
  } else {
    throw std::runtime_error("Неизвестный тип баростата: " + type);
  }
}

int main(int argc, char *argv[]) {
  {
    Timer timer(0);
    timer.start();
    // Проверка аргументов командной строки для имени конфигурационного файла
    if (argc < 2) {
      std::cerr << "Использование: " << argv[0] << " config.json\n";
      return 1;
    }
    std::ifstream configFile(argv[1]);
    if (!configFile) {
      std::cerr << "Не удалось открыть файл конфигурации.\n";
      return 1;
    }
    // Считали конфигурацию
    json config;
    configFile >> config;
    configFile.close();

    // Настройки запуска
    Settings settings(config); // Настройки системы
    cout << settings.getData() << endl;

    // Создаем пул потоков
    ThreadPool threadPool(settings.threads());
    cout << threadPool.getData() << endl;

    // Заводим переменную System, которая хранит все состояние системы
    System sys = System(config, settings);
    cout << sys.getData() << endl;

    // Заводим переменную Backup для бекапов
    BackupManager backupManager = BackupManager(config["backup"], settings);
    cout << backupManager.getData() << endl;
    if (!backupManager.enabled()) {
      // Генерация начальных координат и скоростей
      Generator generator = Generator(settings);
      generator.generateCoords(sys, settings);
      if (settings.structType() != std::string("Custom"))
        generator.generateSpeeds(sys, settings);
    } else {
      if (backupManager.restoreBackup(sys))
        cout << "Восстановление из резервной копии прошло успешно!" << endl;
      else {
        exit(0);
      }
    }

    // Создаём потенциал, термостат и баростат через фабрики:

    // Заводим потенциал (обязателен)
    std::unique_ptr<Potential> potential{nullptr};
    if (!config.contains("potential"))
      throw std::invalid_argument(
          "В файле конфигурации нет настроек используемого потенциала");
    potential = createPotential(config["potential"]);
    cout << potential->getData() << endl;

    // Заводим термостат, если включен
    std::unique_ptr<Thermostat> thermostat{nullptr};
    if (config.contains("thermostat") &&
        config["thermostat"].value("toggle", false)) {
      thermostat =
          createThermostat(config["thermostat"], settings, sys.particleNumber());
    }
    if (thermostat)
      cout << thermostat->getData() << endl;
    else
      cout << "Термостат выключен" << endl;

    // Заводим баростат, если включен
    std::unique_ptr<Barostat> barostat{nullptr};
    if (config.contains("barostat") &&
        config["barostat"].value("toggle", false)) {
      barostat = createBarostat(config["barostat"], settings);
    }
    if (barostat)
      cout << barostat->getData();
    else
      cout << "Баростат выключен" << endl;

    // Заводим переменную Output для вывода в файлы/консоль
    OutputManager outputManager = OutputManager(config["output"], settings);
    cout << outputManager.getData() << endl;

    // Заводим переменную EnsembleManager для управления ансамблями
    std::unique_ptr<EnsembleManager> ensembleManager{nullptr};
    if (config["ensemble"].value("toggle", false)) {
      ensembleManager = std::make_unique<EnsembleManager>(
          config["ensemble"], settings, sys, threadPool, outputManager);
      cout << ensembleManager->getData() << endl;
    }

    CellList cellList;
    // Заводим переменную MDAlgorimths, которая содержит все методы МД
    json macroparams_config = config["macroparams"];
    MDAlgorithms md = MDAlgorithms(settings, sys, backupManager, outputManager,
                                   threadPool,cellList, std::move(ensembleManager),
                                   std::move(thermostat), std::move(barostat), 
                                   std::move(potential), macroparams_config);

    // Записываем конфигурацию запуска
    outputManager.writeModellingProperties(config);

    //================
    // НАЧАЛО РАСЧЕТОВ
    cout << "Начало расчета!" << endl;
    Timer step_timer(0); // Таймер для отслеживания времени выполнения шагов
    step_timer.start();
    md.initialStep(sys); // Начальный (нулевой) шаг
    step_timer.stop();
    cout << sys.getShortData() << endl;
    cout << "Init step done in " << step_timer.elapsed() << endl;

    bool completed = false;
    // Основной цикл МД
    while (!completed && sys.currentStep() < settings.steps()) {
      step_timer.start();
      completed = md.advanceStep(sys);
      step_timer.stop();
      cout << sys.getShortData() << endl;
      cout << "Step " << sys.currentStep() << " - done in "
           << step_timer.elapsed() << endl;
    }
    timer.stop();
    cout << "Time elapsed: " << timer.elapsed() << endl;
    cout << "Завершение расчета" << endl;
  }
  return 0;
}
